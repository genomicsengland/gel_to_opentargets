# Importing Genomics England data into Open Targets

## gel_to_ot.py

Takes an input TSV (probably exported from the `tiering_data` table in LabKey) and generates JSON output (one object per line) which should be valid according to the [opentargets/validator](https://github.com/opentargets/validator).

### Phenotypes

`gel_to_ot.py` needs to have a set of mappings from disease name to EFO term; by default this should be in a file called `phenotypes_text_to_efo.txt`. This can be generated from GEL data as follows:

1. Extract relevant data from `tiering_data` table in LabKey, save as TSV

1. Extract phenotypes as strings from the TSV file

`tail -n +2 tiering_data.tsv | cut -d$'\t' -f 4  | sort -uf > phenotypes.txt`

1. Use opentargets/OnToma to perform the mapping

`ontoma phenotypes.txt phenotypes_text_to_efo.txt`  

1. Strip the OnToma output to create a file (not actually required since `gel_to_ot.py` only uses first 2 fields)

`tail -n +2 phenotypes_text_to_efo.txt | cut -d$'\t' -f 1,2`

## Docker

To run via Docker:

`docker run gel_to_ot --input sample.tsv`

To rebuild the Docker image after changes:

`docker build -t gel_to_ot .`

## Sample data

Note that the data in `sample.tsv` is intended to be representative of the tiering data, but certain fields (e.g. participant_id are placeholders and are not taken from real data)
