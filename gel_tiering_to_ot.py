# Take a TSV file exported from LabKey and convert it to JSON following the OpenTargets "genetic association" schema
# https://github.com/opentargets/json_schema

import argparse
import csv
import json
import logging
import sys
import re
import itertools

SOURCE_ID = "genomics_england_tiering"
PHENOTYPE_MAPPING_FILE = "phenotypes_text_to_efo.txt"
DATABASE_ID = "genomics_england_main_programme"
DATABASE_VERSION = "5.1"  # Change if version changes
SNP_REGEXP = "rs[0-9]{1,}"  # TODO - support more SNP types
GEL_LINK_PREFIX = "http://emb-prod-mre-labkey-01.gel.zone:8080/labkey/query/main-programme/main-programme_v5.1_2018-11-20/executeQuery.view?schemaName=lists&query.queryName=participant&query.participant_id~eq="
FAKE_RS_ID_BASE = 2000000000  # Number to start fake rsIDs at


def main():
    parser = argparse.ArgumentParser(description='Generate Open Targets JSON from an input TSV file')

    parser.add_argument('--input', help="Tiering data TSV input file", required=True, action='store')

    parser.add_argument('--pedigree', help="Pedigree data TSV input file", required=True, action='store')

    parser.add_argument("--log-level", help="Set the log level", action='store', default='WARNING')

    args = parser.parse_args()

    logging.basicConfig()
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.getLevelName(args.log_level))

    logger.info("Reading TSV from " + args.input)

    required_columns = ["sample_id", "phenotype", "db_snp_id", "tier", "genomic_feature_ensembl_id",
                        "genomic_feature_hgnc", "consequence_type", "participant_id", "participant_type",
                        "genotype", "mode_of_inheritance", "rare_diseases_family_id"]

    count = 0

    unknown_phenotypes = set()
    unknown_consequences = set()

    consequence_map = build_consequence_type_to_so_map()
    phenotype_map = read_phenotype_to_efo_mapping(PHENOTYPE_MAPPING_FILE)
    apply_phenotype_mapping_overrides(phenotype_map)

    affected_map = build_affected_map(args.pedigree)

    fake_rs_counter = itertools.count(start=FAKE_RS_ID_BASE)

    with open(args.input) as tsv_file:

        reader = csv.DictReader(tsv_file, delimiter='\t')

        for column in required_columns:
            if column not in reader.fieldnames:
                logger.error("Required column %s does not exist in %s (columns in file are %s)"
                             % (column, args.input, reader.fieldnames))
                sys.exit(1)

        for row in reader:

            # Some rows have several, comma-separated consequence types; only use one in this case
            if ',' in row['consequence_type']:
                row['consequence_type'] = row['consequence_type'].split(',')[0]

            my_instance = build_evidence_strings_object(consequence_map, phenotype_map, affected_map, row,
                                                        fake_rs_counter, unknown_phenotypes, unknown_consequences)
            if my_instance:
                print(json.dumps(my_instance))
                count += 1

    logger.info("Processed %d objects" % count)
    logger.info("Generated synthetic rsIDs for %d entries" % (fake_rs_counter.next() - FAKE_RS_ID_BASE))
    logger.info("%d phenotypes were not found:" % len(unknown_phenotypes))
    for phenotype in unknown_phenotypes:
        logger.info(phenotype)
    logger.info("%d consequences were not found:" % len(unknown_consequences))
    for consequence in unknown_consequences:
        logger.info(consequence)


def build_evidence_strings_object(consequence_map, phenotype_map, affected_map, row, fake_rs_counter,
                                  unknown_phenotypes, unknown_consequences):
    """
    Build a Python object containing the correct structure to match the Open Targets genetics.json schema
    :return:
    """

    logger = logging.getLogger(__name__)

    novel_snp = False
    if not re.match(SNP_REGEXP, row['db_snp_id']):
        original_rsid = row['db_snp_id']
        novel_snp = True
        row['db_snp_id'] = "rs%d" % next(fake_rs_counter)
        logger.info("Record with sample ID %s, Ensembl ID %s and phenotype %s has variant %s which does not match "
                    "the list of allowed types, so generating fake rsID %s" % (
                        row['sample_id'], row['genomic_feature_ensembl_id'],
                        row['phenotype'], original_rsid, row['db_snp_id']))

    logger.debug("Building container object")

    if row['consequence_type'] not in consequence_map:
        unknown_consequences.add(row['consequence_type'])
        return

    functional_consequence = consequence_map[row['consequence_type']]

    phenotype = row['phenotype'].strip()
    if phenotype not in phenotype_map:
        unknown_phenotypes.add(phenotype)
        return

    ontology_term = phenotype_map[phenotype]

    gel_link = GEL_LINK_PREFIX + row['participant_id']

    link_text = build_link_text(row, affected_map, novel_snp)

    score = tier_to_score(row['tier'])

    clinical_significance = tier_to_clinical_significance(row['tier'])

    obj = {
        "sourceID": SOURCE_ID,
        "access_level": "public",
        "validated_against_schema_version": "1.2.8",
        "unique_association_fields": {
            "sample_id": row['sample_id'],
            "participant_id": row['participant_id'],
            "gene": row['genomic_feature_ensembl_id'],
            "phenotype": phenotype,
            "variant": row['db_snp_id']
        },
        "target": {
            "id": "http://identifiers.org/ensembl/" + row['genomic_feature_ensembl_id'],
            "target_type": "http://identifiers.org/cttv.target/gene_evidence",
            "activity": "http://identifiers.org/cttv.activity/loss_of_function"
        },
        "disease": {
            "name": phenotype,
            "id": ontology_term
        },
        "type": "genetic_association",
        "variant": {
            "id": "http://identifiers.org/dbsnp/" + row['db_snp_id'],
            "type": "snp single"
        },
        "evidence": {
            "gene2variant": {
                "is_associated": True,
                "date_asserted": "2018-10-22T23:00:00",
                "provenance_type": {
                    "database": {
                        "id": DATABASE_ID,
                        "version": DATABASE_VERSION
                    }
                },
                "evidence_codes": [
                    "http://identifiers.org/eco/cttv_mapping_pipeline"
                ],
                "urls": [
                    {
                        "nice_name": link_text,
                        "url": gel_link
                    }
                ],
                "functional_consequence": functional_consequence
            },
            "variant2disease": {
                "unique_experiment_reference": "STUDYID_" + row['sample_id'],  # Required by regexp in base.json
                "is_associated": True,
                "date_asserted": "2018-10-22T23:00:00",
                "resource_score": {
                    "type": "probability",
                    "value": score
                },
                "provenance_type": {
                    "database": {
                        "id": DATABASE_ID,
                        "version": DATABASE_VERSION
                    }
                },
                "evidence_codes": [
                    "http://identifiers.org/eco/GWAS"
                ],
                "urls": [
                    {
                        "nice_name": link_text,
                        "url": gel_link
                    }
                ],
                "clinical_significance": clinical_significance
            }
        }
    }

    return obj


def build_consequence_type_to_so_map():
    consequence_to_so_map = {
        "3_prime_UTR_variant": "http://purl.obolibrary.org/obo/SO_0001624",
        "5_prime_UTR_variant": "http://purl.obolibrary.org/obo/SO_0001623",
        "coding_sequence_variant": "http://purl.obolibrary.org/obo/SO_0001580",
        "downstream_gene_variant": "http://purl.obolibrary.org/obo/SO_0001632",
        "feature_elongation": "http://purl.obolibrary.org/obo/SO_0001907",
        "feature_truncation": "http://purl.obolibrary.org/obo/SO_0001906",
        "frameshift_variant": "http://purl.obolibrary.org/obo/SO_0001589",
        "incomplete_terminal_codon_variant": "http://purl.obolibrary.org/obo/SO_0001626",
        "inframe_deletion": "http://purl.obolibrary.org/obo/SO_0001822",
        "inframe_insertion": "http://purl.obolibrary.org/obo/SO_0001821",
        "intergenic_variant": "http://purl.obolibrary.org/obo/SO_0001628",
        "intron_variant": "http://purl.obolibrary.org/obo/SO_0001627",
        "mature_miRNA_variant": "http://purl.obolibrary.org/obo/SO_0001620",
        "missense_variant": "http://purl.obolibrary.org/obo/SO_0001583",
        "NMD_transcript_variant": "http://purl.obolibrary.org/obo/SO_0001621",
        "non_coding_transcript_exon_variant": "http://purl.obolibrary.org/obo/SO_0001792",
        "non_coding_transcript_variant": "http://purl.obolibrary.org/obo/SO_0001619",
        "protein_altering_variant": "http://purl.obolibrary.org/obo/SO_0001818",
        "regulatory_region_ablation": "http://purl.obolibrary.org/obo/SO_0001894",
        "regulatory_region_amplification": "http://purl.obolibrary.org/obo/SO_0001891",
        "regulatory_region_variant": "http://purl.obolibrary.org/obo/SO_0001566",
        "splice_acceptor_variant": "http://purl.obolibrary.org/obo/SO_0001574",
        "splice_donor_variant": "http://purl.obolibrary.org/obo/SO_0001575",
        "splice_region_variant": "http://purl.obolibrary.org/obo/SO_0001630",
        "start_lost": "http://purl.obolibrary.org/obo/SO_0002012",
        "stop_gained": "http://purl.obolibrary.org/obo/SO_0001587",
        "stop_lost": "http://purl.obolibrary.org/obo/SO_0001578",
        "stop_retained_variant": "http://purl.obolibrary.org/obo/SO_0001567",
        "synonymous_variant": "http://purl.obolibrary.org/obo/SO_0001819",
        "TF_binding_site_variant": "http://purl.obolibrary.org/obo/SO_0001782",
        "TFBS_ablation": "http://purl.obolibrary.org/obo/SO_0001895",
        "TFBS_amplification": "http://purl.obolibrary.org/obo/SO_0001892",
        "transcript_ablation": "http://purl.obolibrary.org/obo/SO_0001893",
        "transcript_amplification": "http://purl.obolibrary.org/obo/SO_0001889",
        "upstream_gene_variant": "http://purl.obolibrary.org/obo/SO_0001631"
    }

    return consequence_to_so_map


def read_phenotype_to_efo_mapping(filename):
    phenotype_map = dict()

    with open(filename, 'r') as mapping_file:
        for line in mapping_file:
            line = line.rstrip("\n")
            if line.startswith("#") or line.startswith("query"):
                continue
            parts = line.split("\t")
            phenotype = parts[0].strip()
            ontology_term = parts[1].strip()
            if phenotype and ontology_term:
                phenotype_map[phenotype] = ontology_term

    return phenotype_map


def apply_phenotype_mapping_overrides(phenotype_map):
    phenotype_map["Early onset and familial Parkinson's Disease"] = "http://www.ebi.ac.uk/efo/EFO_0002508"
    phenotype_map["Early onset and familial Parkinson%#27;s Disease"] = "http://www.ebi.ac.uk/efo/EFO_0002508"
    phenotype_map["early onset and familial parkinsons disease"] = "http://www.ebi.ac.uk/efo/EFO_0002508"
    return


def tier_to_score(tier):
    tier = tier.lower()

    if tier == "tier1":
        score = 1
    elif tier == "tier2":
        score = 0.5
    elif tier == "tier3":
        score = 0.1
    else:
        score = 0

    return score


def tier_to_clinical_significance(tier):
    # See https://github.com/opentargets/json_schema/blob/master/opentargets.json for values

    tier = tier.lower()

    if tier == "tier1":
        cs = "Pathogenic"
    elif tier == "tier2":
        cs = "Likely pathogenic"
    elif tier == "tier3":
        cs = "association"
    else:
        cs = "association"

    return cs


def build_affected_map(filename):
    """
    Parse pedigree file to build a map of affected/unaffected status for each participant.
    Returns: Dict of particpant ID, afected status
    """

    affected_map = dict()

    with open(filename) as pedigree_file:

        reader = csv.DictReader(pedigree_file, delimiter='\t')

        for row in reader:

            if not row['member_participant_id']:
                continue

            if not row['affection_status']:
                continue

            affected_map[row['member_participant_id']] = row['affection_status'].lower()

    return affected_map


def build_link_text(row, affected_map, novel_snp):
    """
    Build text that is displayed on participant link, e.g.
    GeL TIER2  monoallelic_not_imprinted variant for family G105275; participant 113001766 (affected) is heterozygous.
    :param row: dict of columns in current evidence line
    :param affected_map: dict of participant_id:affected status
    :param novel_snp: whether this SNP is novel or not
    :return: String of text
    """
    pid = row['participant_id']

    affected = affected_map.get(pid, "unknown")

    text = "GeL {tier} {mode} {snp_status} variant for family {family}; participant {participant} ({affected}) is {genotype}".format(
        tier=row['tier'],
        genotype=row['genotype'],
        snp_status="novel" if novel_snp else "known",
        family=row['rare_diseases_family_id'],
        participant=pid,
        affected=affected,
        mode=row['mode_of_inheritance'])

    return text


if __name__ == '__main__':
    sys.exit(main())
