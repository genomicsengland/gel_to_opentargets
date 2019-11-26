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
DATABASE_VERSION = "8"  # Change if version changes
ASSERTION_DATE = "2019-11-28T23:00:00"  # Change to date of data release
LABKEY_QUESTIONNAIRE_LINK_TEMPLATE = "http://emb-prod-mre-labkey-01.gel.zone:8080/labkey/query/main-programme/main-programme_v8_2019-11-28/executeQuery.view?schemaName=lists&query.queryName=gmc_exit_questionnaire&query.participant_id~eq={participant}"
SCHEMA_VERSION = "1.6.0"  # Open Targets JSON schema version


def main():
    parser = argparse.ArgumentParser(description='Generate Open Targets exit questionnaire JSON from an input TSV file')

    parser.add_argument('--input', help="Questionnaire data TSV input file", required=True, action='store')

    parser.add_argument('--hgnc_to_ensembl', help="File containing a list of HGNC symbol to Ensembl gene ID mappings",
                        required=True, action='store')

    # TODO - add argument for rare_diseases_participant_disease file

    parser.add_argument('--filter_participants', help="List of participants to filter out", required=False, action='store')

    parser.add_argument("--log-level", help="Set the log level", action='store', default='WARNING')

    args = parser.parse_args()

    logging.basicConfig()
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.getLevelName(args.log_level))

    required_columns = ["participant_id", "family_id", "chromosome", "position", "reference",
                        "alternate", "acmg_classification", "actionability", "case_solved_family", "gene_name"]

    count = 0

    unknown_phenotypes = set()

    phenotype_map = gel_utils.read_phenotype_to_efo_mapping(PHENOTYPE_MAPPING_FILE)
    gel_utils.apply_phenotype_mapping_overrides(phenotype_map)

    # Read tiering data to get variant:ensembl gene mapping
    hgnc_to_ensembl = read_hgnc_to_ensembl_mapping(args.hgnc_to_ensembl)
    logger.debug('Read {} HGNC:Ensembl mappings'.format(len(hgnc_to_ensembl)))

    if args.filter_participants:
        participants_to_filter = gel_utils.read_participants_to_filter(args.filter_participants, logger)
    else:
        participants_to_filter = list()

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
            if row['participant_id'] in participants_to_filter:
                continue

            my_instance = build_evidence_strings_object(row, phenotype_map, unknown_phenotypes, hgnc_to_ensembl)

            if my_instance:
                print(json.dumps(my_instance))
                count += 1

    logger.info("Processed {} objects".format(count))
    logger.info("{} phenotypes were not found:".format(len(unknown_phenotypes)))
    for phenotype in unknown_phenotypes:
        logger.info(phenotype)

def build_evidence_strings_object(row, phenotype_map, unknown_phenotypes, hgnc_to_ensembl):
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

    score = 1  # TODO different score based on positive or negative result - e.g. 0 or skip entirely if phenotypes_solved is "no"?

    clinical_significance = row['acmg_classification']

    if row['gene_name'] == 'NA':
        # TODO record number of NAs / missed lookups
        return

    # Only use first gene name if there are multiples separated by ;
    gene_name = row['gene_name'].split(';')[0]

    if gene_name not in hgnc_to_ensembl:
        print "No Ensembl ID found for HGNC symbol " + row['gene_name']
        # TODO what to do here?
    else:
        gene = hgnc_to_ensembl[gene_name]

    # Build composite variant
    variant = ':'.join((row['chromosome'], row['position'], row['reference'], row['alternate']))  # matches format in map

    # Link to LabKey based on participant only
    gel_link = LABKEY_QUESTIONNAIRE_LINK_TEMPLATE.format(participant=participant_id)

    link_text = build_link_text(row)

    obj = {
        "sourceID": SOURCE_ID,
        "access_level": "public",
        "validated_against_schema_version": SCHEMA_VERSION,
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
            "id": variant,
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


def read_hgnc_to_ensembl_mapping(hgnc_to_ensembl_file_name):
    """
    Build a map of HGNC symbols (used in GEL questionnaire data) to Ensembl gene IDs (used in Open Targets).
    :param hgnc_to_ensembl_file_name: Name of mapping file.
    :return: Map of HGNC to Ensembl identifiers.
    """
    hgnc_to_ensembl = {}

    with open(hgnc_to_ensembl_file_name, 'r') as mapping_file:
        for line in mapping_file:
            (hgnc, ensembl) = line.split()
            hgnc_to_ensembl[hgnc] = ensembl

    return hgnc_to_ensembl


if __name__ == '__main__':
    sys.exit(main())
