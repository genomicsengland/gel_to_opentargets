# Importing Genomics England data into Open Targets

## gel_to_ot.py

Takes an input TSV (probably exported from the `tiering_data` table in LabKey) and generates JSON output (one object per line) which should be valid according to the [opentargets/validator](https://github.com/opentargets/validator).

### Phenotypes

`gel_to_ot.py` needs to have a set of mappings from disease name to EFO term; by default this should be in a file called `phenotypes_text_to_efo.txt`. This can be generated from GEL data as follows:

1. Extract relevant data from `tiering_data` table in LabKey, save as TSV

2. Extract phenotypes as strings from the TSV file

`tail -n +2 tiering_data.tsv | cut -d$'\t' -f 4  | sort -uf > phenotypes.txt`

3. Use opentargets/OnToma to perform the mapping

`ontoma phenotypes.txt phenotypes_text_to_efo.txt`  

4. Strip the OnToma output to create a file (not actually required since `gel_to_ot.py` only uses first 2 fields)

`tail -n +2 phenotypes_text_to_efo.txt | cut -d$'\t' -f 1,2`

## Pedigree

Extract the `rare_diseases_pedigree_member` table from LabKey and save as a TSV. This is used to determine affected/unaffected status.

### Novel variants

Currently the Open Targets schema requires that all genetic associations have variants with a dbSNP rsID assigned. This is not possible in the case of novel variants. As a workaround, placeholder rsIDs, starting at rs2000000000 are assigned to novel variants during the import process.

## Usage

`python gel_tiering_to_ot.py --input sample.tsv --pedigree sample_pedigree.tsv`

## Docker

To run via Docker:

`docker run gel_to_ot --input test_data.tsv --pedigree sample_pedigree.tsv`

To rebuild the Docker image after changes:

`docker build -t gel_to_ot .`

## Test data

Note that the data in `test_data.tsv` is intended to be representative of the tiering data, but certain fields (e.g. participant_id are placeholders and are not taken from real data)
