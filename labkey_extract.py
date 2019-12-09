# Extract data required for Open Targets from LabKey

# See https://cnfl.extge.co.uk/display/GERE/Python+LabKey+API

# Note need to populate .netrc and do
#   module load python/2.7.12

import labkey
import pandas as pd

LABKEY_SERVER = "emb-prod-mre-labkey-01.gel.zone:8080"
RELEASE_PATH = "main-programme/main-programme_v8_2019-11-28"

API_KEY = "" # TODO fill in when available

def labkey_to_api(table, filename, filter_array, server_context):

    print "Extracting " + table + " to " + filename
    
    question_results = labkey.query.select_rows(server_context = server_context,
                                                schema_name = 'lists',
                                                query_name = table,
                                                filter_array = filter_array,
                                                max_rows = 1000000)

    df = pd.DataFrame(question_results["rows"])

    df.to_csv(filename, sep='\t', encoding='utf-8')



print "Connecting to LabKey on " + LABKEY_SERVER + " " + RELEASE_PATH

server_context = labkey.utils.create_server_context(LABKEY_SERVER, RELEASE_PATH, 'labkey', use_ssl=False, api_key=API_KEY)

labkey_to_api("gmc_exit_questionnaire", "gmc_exit_questionnaire.tsv", [], server_context)

labkey_to_api("rare_diseases_pedigree_member", "rare_diseases_pedigree_member.tsv", [], server_context)

filter_array = [ labkey.query.QueryFilter('tier', 'TIER1;TIER2', 'in') ]

labkey_to_api("tiering_data", "tiering_data.tsv", filter_array, server_context)

labkey_to_api("rare_diseases_participant_disease", "rare_diseases_participant_disease.tsv", [], server_context)

