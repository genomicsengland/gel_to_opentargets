"""
Utility functions for converting GEL data to Open Targets.
"""


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
    """
    Modify a map of phenotype names to EFO URLs to add new mappings or override existing ones.
    :param phenotype_map: Dict of phenotype names to EFO URLs
    :return: Nothing (the dict is passed by value and modified in place)
    """
    phenotype_map["Early onset and familial Parkinson's Disease"] = "http://www.ebi.ac.uk/efo/EFO_0002508"
    phenotype_map["Early onset and familial Parkinson%#27;s Disease"] = "http://www.ebi.ac.uk/efo/EFO_0002508"
    phenotype_map["early onset and familial parkinsons disease"] = "http://www.ebi.ac.uk/efo/EFO_0002508"
    return


def build_consequence_type_to_so_map():
    """
    Build a cache of consequence types (e.g. "stop_lost") to OBO URLS (e.g. "http://purl.obolibrary.org/obo/SO_0002012")
    :return: Dict containing URLs for each consequence type
    """
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
