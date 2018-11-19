# Take a TSV file exported from LabKey and convert it to JSON following the OpenTargets "genetic association" schema
# https://github.com/opentargets/json_schema

import argparse
import csv
import json
import logging
import sys
import re

SOURCE_ID = "genomics_england_tiering"
PHENOTYPE_MAPPING_FILE = "phenotypes_text_to_efo.txt"
DATABASE_ID = "genomics_england_tiering"  # TODO Change to GEL Main Programme when working
DATABASE_VERSION = "1.0"  # Change if version changes
SNP_REGEXP = "rs[0-9]{1,}"  # TODO - support more SNP types
GEL_LINK_PREFIX = "https://opencga-embassy.gel.zone/iva/#browser/reopencga@100k_genomes_grch37_germline/RD37?gene="


def main():
    parser = argparse.ArgumentParser(description='Generate Open Targets JSON from an input TSV file')

    parser.add_argument('--input', help="TSV input file", required=True, action='store')

    parser.add_argument("--log-level", help="Set the log level", action='store', default='WARNING')

    args = parser.parse_args()

    logging.basicConfig()
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.getLevelName(args.log_level))

    logger.info("Reading TSV from " + args.input)

    required_columns = ["sample_id", "phenotype", "db_snp_id", "tier", "genomic_feature_ensembl_id",
                        "genomic_feature_hgnc", "consequence_type"]

    count = 0

    consequence_map = build_consequence_type_to_so_map()
    phenotype_map = read_phenotype_to_efo_mapping(PHENOTYPE_MAPPING_FILE)

    with open(args.input) as tsv_file:

        reader = csv.DictReader(tsv_file, delimiter='\t')

        for column in required_columns:
            if column not in reader.fieldnames:
                logger.error("Required column %s does not exist in %s (columns in file are %s)"
                             % (column, args.input, reader.fieldnames))
                sys.exit(1)

        for row in reader:

            # Some rows have several, comma-separated consequence types; split these out into separate rows
            for consequence_type in row['consequence_type'].split(","):
                row['consequence_type'] = consequence_type
                my_instance = build_evidence_strings_object(consequence_map, phenotype_map, row)
                if my_instance:
                    print(json.dumps(my_instance))
                    count += 1

    logger.info("Processed %d objects" % count)


def build_evidence_strings_object(consequence_map, phenotype_map, row):
    """
    Build a Python object containing the correct structure to match the Open Targets genetics.json schema
    :return:
    """

    logger = logging.getLogger(__name__)

    if not re.match(SNP_REGEXP, row['db_snp_id']):
        logger.info("Record with sample ID %s, Ensembl ID %s and phenotype %s has variant %s which does not match "
                    "the list of allowed types, ignoring" % (row['sample_id'], row['genomic_feature_ensembl_id'],
                                                             row['phenotype'], row['db_snp_id']))
        return

    logger.debug("Building container object")

    if row['consequence_type'] not in consequence_map:
        logger.error("Can't find consequence type mapping for %s, skipping this row" % row['consequence_type'])
        return

    functional_consequence = consequence_map[row['consequence_type']]

    phenotype = row['phenotype'].strip()
    if phenotype not in phenotype_map:
        logger.error("Can't find phenotype mapping for %s, skipping this row" % phenotype)
        return

    ontology_term = phenotype_map[phenotype]

    gel_link = GEL_LINK_PREFIX + row['db_snp_id']

    score = tier_to_score(row['tier'])

    obj = {
        "sourceID": SOURCE_ID,
        "access_level": "public",
        "validated_against_schema_version": "1.2.8",
        "unique_association_fields": {
            "sample_id": row['sample_id'],
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
                        "nice_name": "Sourced from Genomics England tiering data",
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
                        "nice_name": "Sourced from Genomics England tiering data",
                        "url": "https://www.genomicsengland.co.uk/"  # TODO change for more specific URL
                    }
                ]
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


if __name__ == '__main__':
    sys.exit(main())
