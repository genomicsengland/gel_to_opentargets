# Take a TSV file exported from LabKey and convert it to JSON following the OpenTargets "genetic association" schema
# https://github.com/opentargets/json_schema

import argparse
import csv
import json
import logging
import sys
import gel_utils

SOURCE_ID = "genomics_england_tiering"
PHENOTYPE_MAPPING_FILE = "phenotypes_text_to_efo.txt"
DATABASE_ID = "genomics_england_main_programme"
DATABASE_VERSION = "8"  # Change if version changes
ASSERTION_DATE = "2019-11-28T23:00:00"  # Change to date of data release
LABKEY_PARTICIPANT_LINK_TEMPLATE = "http://emb-prod-mre-labkey-01.gel.zone:8080/labkey/query/main-programme/main-programme_v8_2019-11-28/executeQuery.view?schemaName=lists&query.queryName=participant&query.participant_id~eq={participant_id}"
SCHEMA_VERSION = "1.6.0"  # Open Targets JSON schema version


def main():
    parser = argparse.ArgumentParser(description='Generate Open Targets JSON from an input TSV file')

    parser.add_argument('--input', help="Tiering data TSV input file", required=True, action='store')

    parser.add_argument('--pedigree', help="Pedigree data TSV input file", required=True, action='store')

    parser.add_argument('--filter_participants', help="List of participants to filter out", required=False, action='store')

    parser.add_argument("--log-level", help="Set the log level", action='store', default='WARNING')

    args = parser.parse_args()

    logging.basicConfig()
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.getLevelName(args.log_level))

    logger.info("Reading TSV from " + args.input)

    required_columns = ["sample_id", "phenotype", "db_snp_id", "tier", "ensembl_id",
                        "genomic_feature_hgnc", "consequence_type", "participant_id", "participant_type",
                        "genotype", "mode_of_inheritance", "rare_diseases_family_id"]

    count = 0

    unknown_phenotypes = set()
    unknown_consequences = set()

    consequence_map = gel_utils.build_consequence_type_to_so_map()
    phenotype_map = gel_utils.read_phenotype_to_efo_mapping(PHENOTYPE_MAPPING_FILE)
    gel_utils.apply_phenotype_mapping_overrides(phenotype_map)

    affected_map = build_affected_map(args.pedigree)

    if args.filter_participants:
        participants_to_filter = gel_utils.read_participants_to_filter(args.filter_participants, logger)
    else:
        participants_to_filter = list()

    with open(args.input) as tsv_file:

        reader = csv.DictReader(tsv_file, delimiter='\t')

        for column in required_columns:
            if column not in reader.fieldnames:
                logger.error("Required column %s does not exist in %s (columns in file are %s)"
                             % (column, args.input, reader.fieldnames))
                sys.exit(1)

        for row in reader:

            if row['participant_id'] in participants_to_filter:
                continue

            # Some rows have several, comma-separated consequence types; only use one in this case
            if ',' in row['consequence_type']:
                row['consequence_type'] = row['consequence_type'].split(',')[0]

            my_instance = build_evidence_strings_object(consequence_map, phenotype_map, affected_map, row,
                                                        unknown_phenotypes, unknown_consequences)
            if my_instance:
                print(json.dumps(my_instance))
                count += 1

    logger.info("Processed %d objects" % count)
    logger.info("%d phenotypes were not found:" % len(unknown_phenotypes))
    for phenotype in unknown_phenotypes:
        logger.info(phenotype)
    logger.info("%d consequences were not found:" % len(unknown_consequences))
    for consequence in unknown_consequences:
        logger.info(consequence)


def build_evidence_strings_object(consequence_map, phenotype_map, affected_map, row,
                                  unknown_phenotypes, unknown_consequences):
    """
    Build a Python object containing the correct structure to match the Open Targets genetics.json schema
    :return:
    """

    logger = logging.getLogger(__name__)

    if row['consequence_type'] not in consequence_map:
        unknown_consequences.add(row['consequence_type'])
        return

    functional_consequence = consequence_map[row['consequence_type']]

    phenotype = row['phenotype'].strip()
    if phenotype not in phenotype_map:
        unknown_phenotypes.add(phenotype)
        return

    ontology_term = phenotype_map[phenotype]

    gel_link = LABKEY_PARTICIPANT_LINK_TEMPLATE.format(participant_id=row['participant_id'])

    link_text = build_link_text(row, affected_map)

    score = tier_to_score(row['tier'])

    clinical_significance = tier_to_clinical_significance(row['tier'])

    obj = {
        "sourceID": SOURCE_ID,
        "access_level": "public",
        "validated_against_schema_version": SCHEMA_VERSION,
        "unique_association_fields": {
            "sample_id": row['sample_id'],
            "participant_id": row['participant_id'],
            "gene": row['ensembl_id'],
            "phenotype": phenotype,
            "variant": row['db_snp_id']
        },
        "target": {
            "id": "http://identifiers.org/ensembl/" + row['ensembl_id'],
            "target_type": "http://identifiers.org/cttv.target/gene_evidence",
            "activity": "http://identifiers.org/cttv.activity/loss_of_function"
        },
        "disease": {
            "name": phenotype,
            "id": ontology_term
        },
        "type": "genetic_association",
        "variant": {
            "id": "http://identifiers.org/dbsnp/rs0",
            "type": "snp single"
        },
        "evidence": {
            "gene2variant": {
                "is_associated": True,
                "date_asserted": ASSERTION_DATE,
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
                "date_asserted": ASSERTION_DATE,
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

    # for the tiering data we can only infer an association, not whether it's pathogenic or not
    return "association"


def build_affected_map(filename):
    """
    Parse pedigree file to build a map of affected/unaffected status for each participant.
    Returns: Dict of participant ID, affected status
    """

    affected_map = dict()

    with open(filename) as pedigree_file:

        reader = csv.DictReader(pedigree_file, delimiter='\t')

        for row in reader:

            if not row['participant_id']:
                continue

            if not row['affection_status']:
                continue

            affected_map[row['participant_id']] = row['affection_status'].lower()

    return affected_map


def build_link_text(row, affected_map):
    """
    Build text that is displayed on participant link, e.g.
    GeL TIER2  monoallelic_not_imprinted variant for family G105275; participant 113001766 (affected) is heterozygous.
    :param row: dict of columns in current evidence line
    :param affected_map: dict of participant_id:affected status
    :return: String of text
    """
    pid = row['participant_id']

    affected = affected_map.get(pid, "unknown")

    text = "GeL {tier} {mode} variant for family {family}; participant {participant} ({affected}) is {genotype}".format(
        tier=row['tier'],
        genotype=row['genotype'].replace('_', ' '),
        family=row['rare_diseases_family_id'],
        participant=pid,
        affected=affected,
        mode=row['mode_of_inheritance'].replace('_', ' '))

    return text


if __name__ == '__main__':
    sys.exit(main())
