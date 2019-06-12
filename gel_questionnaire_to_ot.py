# Take a TSV file of exit questionnaire data exported from LabKey and convert it to JSON following
# the OpenTargets "genetic association" schema https://github.com/opentargets/json_schema

import argparse
import csv
import json
import logging
import sys
import gel_utils

SOURCE_ID = "genomics_england_questionnaire"
PHENOTYPE_MAPPING_FILE = "phenotypes_text_to_efo.txt"
DATABASE_ID = "genomics_england_main_programme"
DATABASE_VERSION = "6"  # Change if version changes
ASSERTION_DATE = "2019-02-28T23:00:00"  # Change to date of data release
LABKEY_QUESTIONNAIRE_LINK_TEMPLATE = "http://emb-prod-mre-labkey-01.gel.zone:8080/labkey/query/main-programme/main-programme_v6_2019-02-28/executeQuery.view?schemaName=lists&query.queryName=gmc_exit_questionnaire&query.variant_details~eq={variant}&query.participant_id~eq={participant}"


def main():
    parser = argparse.ArgumentParser(description='Generate Open Targets exit questionnaire JSON from an input TSV file')

    parser.add_argument('--input', help="Questionnaire data TSV input file", required=True, action='store')

    parser.add_argument('--tiering', help="Tiering data TSV input file, required for variant:gene mapping",
                        required=True, action='store')

    parser.add_argument("--log-level", help="Set the log level", action='store', default='WARNING')

    args = parser.parse_args()

    logging.basicConfig()
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.getLevelName(args.log_level))

    required_columns = ["participant_id", "family_id", "phenotypes_explained", "variant_details",
                        "acmg_classification", "actionability", "case_solved_family"]

    count = 0

    unknown_phenotypes = set()

    phenotype_map = gel_utils.read_phenotype_to_efo_mapping(PHENOTYPE_MAPPING_FILE)
    gel_utils.apply_phenotype_mapping_overrides(phenotype_map)

    # Read tiering data to get variant:ensembl gene mapping
    variant_to_gene = read_variant_to_gene_map_from_tiering(args.tiering)
    logger.debug('Read {} variant:gene mappings'.format(len(variant_to_gene)))

    unknown_variants = set()

    logger.info("Reading TSV from " + args.input)

    with open(args.input) as question_tsv_file:

        reader = csv.DictReader(question_tsv_file, delimiter='\t')

        for column in required_columns:
            if column not in reader.fieldnames:
                logger.error(
                    "Required column {} does not exist in {} (columns in file are {})".format(column, args.input,
                                                                                              reader.fieldnames))
                sys.exit(1)

        for row in reader:

            my_instance = build_evidence_strings_object(row, phenotype_map, unknown_phenotypes, variant_to_gene,
                                                        unknown_variants)
            if my_instance:
                print(json.dumps(my_instance))
                count += 1

    logger.info("Processed {} objects".format(count))
    logger.info("{} phenotypes were not found:".format(len(unknown_phenotypes)))
    for phenotype in unknown_phenotypes:
        logger.info(phenotype)
    logger.info("{} variants could not be mapped to genes and were skipped:".format(len(unknown_variants)))


def build_evidence_strings_object(row, phenotype_map, unknown_phenotypes, variant_to_gene, unknown_variants):
    """
    Build a Python object containing the correct structure to match the Open Targets genetics.json schema
    :return:
    """

    logger = logging.getLogger(__name__)

    participant_id = row['participant_id']

    phenotype = row['phenotypes_explained'].strip()
    if phenotype not in phenotype_map:
        unknown_phenotypes.add(phenotype)
        return

    ontology_term = phenotype_map[phenotype]

    score = 1  # TODO different score based on positibve or negative result - e.g. 0 or skip entirely if phenotypes_solved is "no"?

    clinical_significance = row['acmg_classification']  # TODO check if it's in the allowed enum

    variant = row['variant_details']
    if variant not in variant_to_gene:
        logger.error("No gene found for variant {}".format(variant))
        unknown_variants.add()
        return
    else:
        gene = variant_to_gene[variant]

    gel_link = LABKEY_QUESTIONNAIRE_LINK_TEMPLATE.format(variant=variant, participant=participant_id)

    link_text = build_link_text(row)

    # TODO think about unique association fields

    obj = {
        "sourceID": SOURCE_ID,
        "access_level": "public",
        "validated_against_schema_version": "1.6.0",
        "unique_association_fields": {
            "participant_id": participant_id,
            "gene": gene,
            "phenotype": phenotype,
            "variant": variant
        },
        "target": {
            "id": "http://identifiers.org/ensembl/" + gene,
            "target_type": "http://identifiers.org/cttv.target/gene_evidence",
            "activity": "http://identifiers.org/cttv.activity/loss_of_function"
        },
        "disease": {
            "name": phenotype,
            "id": ontology_term
        },
        "type": "genetic_association",
        "variant": {
            "id": variant,  # TODO wil this validate?
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
                ]
                # functional consequence goes here if needed
            },
            "variant2disease": {
                # TODO check that this is actually unique
                "unique_experiment_reference": participant_id + variant + phenotype,
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


def build_link_text(row):
    """
    Build text that is displayed on participant link, e.g.
    GeL variant for family 1234, participant 4567 case solved, actionable (pathogenic variant)
    :return: String of text
    """

    case_solved = "case solved" if row['case_solved_family'] == 'yes' else 'case not solved'
    actionable = "actionable" if row['actionability'] == 'yes' else 'not actionable'
    classification = row['acmg_classification'].replace('_', ' ')

    text = "GeL variant for family {family}; participant {participant} {solved} {actionable} ({classification})".format(
        family = row['family_id'],
        participant = row['participant_id'],
        solved = case_solved,
        actionable = actionable,
        classification = classification)

    return text


def read_variant_to_gene_map_from_tiering(tiering_file_name):
    variant_to_gene = {}

    with open(tiering_file_name) as tiering_tsv_file:
        reader = csv.DictReader(tiering_tsv_file, delimiter='\t')

        for row in reader:
            variant = ':'.join((row['chromosome'], row['position'], row['reference'], row['alternate']))
            variant_to_gene[variant] = row['ensembl_id']

    return variant_to_gene


if __name__ == '__main__':
    sys.exit(main())
