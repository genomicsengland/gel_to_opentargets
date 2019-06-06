
# TODO document functions

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
