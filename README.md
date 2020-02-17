# Importing Genomics England data into Open Targets

## Data sources

Two Genomics England data sources are currently supported:
 * Rare disease tiering data from the `tiering_data` table in LabKey
 * Clinical exit questionnaire data from the `gmc_exit_questionnaire` table in LabKey
 
This repo contains scripts which will convert the TSV file into JSON (one object per line) according to the Open Targets [JSON schema](https://github.com/opentargets/json_schema) which should be valid according to the [opentargets/validator](https://github.com/opentargets/validator)

There are separate scripts for converting the tiering data and the questionnaire data.

## Preparation

### Phenotypes

Both scripts need to have a set of mappings from disease name to EFO term; by default this should be in a file called `phenotypes_text_to_efo.txt`. _A copy of this is supplied in this repo_, but if you need to generate a new one, instructions are in the "Generating phenotype to EFO mappings" section below.

### Pedigree

Extract the `rare_diseases_pedigree_member` table from LabKey and save as a TSV. This is used to determine affected/unaffected status for the tiering data.

### Participant to disease mapping

Extract the `rare_diseases_participant_disease` table from LabKey and save as a TSV. This is used to assign phenotypes to participants in the questionnaire data import process.

### Filtering participants

If some participants need to be removed (e.g. the TracerX set), create a file containing them, one per line, and specify this via the --filter-participants to both `gel_tiering_to_ot.py` and `gel_questionnaire_to_ot.py`. 
A TSV file containing the TracerX participants can be created in LabKey from the `normalised_consent` column in the `participant` table via http://emb-prod-mre-labkey-01.gel.zone:8080/labkey/query/main-programme/main-programme_v6_2019-02-28/executeQuery.view?schemaName=lists&query.queryName=participant&query.columns=participant_id%2Cnormalised_consent_form&query.normalised_consent_form~eq=cohort-tracerx - note the prefix will change per data release.

### rsIDs

Currently the Open Targets schema requires that all genetic associations have variants with a dbSNP rsID assigned. This is no longer possible (as of release 8) as the IDs are not supplied by GEL. A placeholder ID of `rs0` is assigned in this case.

## Usage

### Extracting data from LabKey

This can be done in 2 ways, either automatically using a Python script:

`python labkey_extract.py`

or manually - see [Manually extracting data from LabKey](#manually-extracting-data-from-labkey).

### Rare disease tiering data

`python gel_tiering_to_ot.py --input test_data.tsv --pedigree sample_pedigree.tsv > tiering.json`

### GMC Exit Questionnaire data

`python gel_questionnaire_to_ot.py --input gmc_exit_questionnaire.tsv python gel_questionnaire_to_ot.py --input questionnaire_test_data.tsv --hgnc_to_ensembl hgnc_to_ensembl.txt --disease_file sample_rare_diseases_participant_disease.tsv > questionnaire.json`

## Docker

A `Dockerfile` is provided which will build an image from which both scripts can be run.

To build the Docker image (or rebuild it after changes):

`docker build -t gel_to_ot .`

To run via Docker:

`docker run gel_to_ot --input test_data.tsv --pedigree sample_pedigree.tsv`


## Test data

Note that the data in `test_data.tsv` is intended to be representative of the tiering data, but certain fields (e.g. participant_id are placeholders and are not taken from real data). Similarly `questionnaire_test_data.tsv` is representative but synthetic.

### Generating phenotype to EFO mappings

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

## Manually extracting data from LabKey

Log in to LabKey and navigate to the appropriate Main Programme release. 

Extract each of the following tables to a TSV file. All data in the table should be used with no filters. Navigate to the appropriate table and click on the "Export" button. Then select the "Text" tab, check that it's set to "Tab separated" and then click "Export to text".
 * `gmc_exit_questionnaire`
 * `rare_diseases_pedigree_member`
 * `rare_diseases_participant_disease`

The final table is `tiering_data`, however we only require rows with the Tier column set to `TIER1` or `TIER2`. This filter can be created as follows: navigate to the `tiering_data` table, click on the "Tier" column heading, select "Filter". A filter window appears. On the "Choose values" tab, select _only_ `TIER1` and `TIER2`. Click "OK" and the filter will be applied. You can check this by looking at the number of rows displayed in the top right-hand corner. It should go from approximately 24 million (unfiltered) to 275,000 (filtered).
Export the filtered table to TSV as described above.

Note that the following URL can be used to access the filtered tiering data directly: http://emb-prod-mre-labkey-01.gel.zone:8080/labkey/query/main-programme/main-programme_v8_2019-11-28/executeQuery.view?schemaName=lists&query.queryName=tiering_data&query.tier~neqornull=TIER3 (the URL may need to be modifed for the release of choice)

You should end up with 4 TSV files (`gmc_exit_questionnaire.tsv`, `rare_diseases_pedigree_member.tsv`, `rare_diseases_participant_disease.tsv`, `tiering_data.tsv`).


