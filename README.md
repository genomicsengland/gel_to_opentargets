# Importing Genomics England data into Open Targets

## Data sources

Two Genomics England data sources are currently supported:
 * Rare disease tiering data from the `tiering_data` table in LabKey
 * Clinical exit questionnaire data from the `gmc_exit_questionnaire` table in LabKey
 
This repo contains scripts which will convert the TSV file into JSON (one object per line) according to the Open Targets [JSON schema](https://github.com/opentargets/json_schema) which should be valid according to the [opentargets/validator](https://github.com/opentargets/validator)

There are separate scripts for converting the tiering data and the questionnaire data.

## Preparation

### Phenotypes

Both scripts need to have a set of mappings from disease name to EFO term; by default this should be in a file called `phenotypes_text_to_efo.txt`. This can be generated from GEL data as follows:

1. Extract relevant data from `tiering_data` table in LabKey, save as TSV. Only rows from `tiering_data` with the `Tier` field having a value of `TIER1` or `TIER2` should be used.

2. Extract relevant data from `rare_diseases_participant_disease` table in LabKey, save as TSV. Currently the whole table is used.

3. Extract phenotypes as strings from the tiering TSV file file (including removing blank lines)

`tail -n +2 tiering_data.tsv | cut -d$'\t' -f 6 > tiering_phenotypes.txt`

4. Extract phenotypes from the exit questionnaire data

`tail -n +2 rare_diseases_participant_disease.tsv | cut -d$'\t' -f 8 > question_phenotypes.txt`

5. Produce a unique set of phenotypes and remove blank lines

`sort -uf tiering_phenotypes.txt question_phenotypes.txt | sed '/^$/d' > phenotypes.txt`

6. Use [opentargets/OnToma](https://github.com/opentargets/OnToma) to perform the mapping

`ontoma phenotypes.txt phenotypes_text_to_efo.txt`  

7. Strip the OnToma output to create a file (not actually required since `gel_to_ot.py` only uses first 2 fields)

`tail -n +2 phenotypes_text_to_efo.txt | cut -d$'\t' -f 1,2`

8. Note that there are some manual phenotype mappings in `gel_utils.py` in the `apply_phenotype_mapping_overrides` function.


### Pedigree

Extract the `rare_diseases_pedigree_member` table from LabKey and save as a TSV. This is used to determine affected/unaffected status for the tiering data.

### Filtering participants

If some participants need to be removed (e.g. the TracerX set), create a file containing them, one per line, and specify this via the --filter-participants to both `gel_tiering_to_ot.py` and `gel_questionnaire_to_ot.py`. 
A TSV file containing the TracerX participants can be created in LabKey from the `normalised_consent` column in the `particiapant` table via http://emb-prod-mre-labkey-01.gel.zone:8080/labkey/query/main-programme/main-programme_v6_2019-02-28/executeQuery.view?schemaName=lists&query.queryName=participant&query.columns=participant_id%2Cnormalised_consent_form&query.normalised_consent_form~eq=cohort-tracerx - note the prefix will change per data release.

### Novel variants

Currently the Open Targets schema requires that all genetic associations have variants with a dbSNP rsID assigned. This is not possible in the case of novel variants. As a workaround, placeholder rsIDs, starting at rs2000000000 are assigned to novel variants during the import process.

## Usage

### Rare disease tiering data

`python gel_tiering_to_ot.py --input sample.tsv --pedigree sample_pedigree.tsv`

### GMC Exit Questionnaire data

`python gel_questionnaire_to_ot.py --input gmc_exit_questionnaire_v6.tsv --tiering tiering_data_v6_t1t2only.tsv`

Note that the tiering data is required for the exit questionnaire conversion; Open Targets requires variants to be mapped to genes, and the questionnaire table does not include Ensembl gene identifiers. The tiering data does, so this is used to generate a lookup table which associates variant locations with genes.

## Docker

A `Dockerfile` is provided which will build an image from which both scripts can be run.

To build the Docker image (or rebuild it after changes):

`docker build -t gel_to_ot .`

To run via Docker:

`docker run gel_to_ot --input test_data.tsv --pedigree sample_pedigree.tsv`


## Test data

Note that the data in `test_data.tsv` is intended to be representative of the tiering data, but certain fields (e.g. participant_id are placeholders and are not taken from real data). Similarly `questionnaire_test_data.tsv` is representative but synthetic.
