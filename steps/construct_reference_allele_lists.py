from typing import Dict, List, Union
import json


def construct_reference_allele_lists(config:Dict, **kwargs) -> None:
    """
    This function takes a dataset and generate an allele list and associated sequence lists for all Class I loci contained within it.

    Args:
        locus (str): the locus to be parsed
        verbose (bool): whether specific information is output to the terminal, for large sequence sets this can be overwhelming and significantly slow down the function
    Returns:

    """
    
    locus = kwargs['locus']

    input_filename = f"output/processed_data/protein_alleles/{locus.lower()}.json"

    with open(input_filename, 'r') as alleles_file:
        alleles = json.load(alleles_file)

    reference_alleles = {
        'allele_groups': {},
        'reference_alleles': [],
    }

    allele_groups = {}




    for allele in alleles:
        allele_group = '_'.join(allele.split('_')[:3]) 
        if allele_group not in reference_alleles['allele_groups']:
            reference_alleles['allele_groups'][allele_group] = allele
            reference_alleles['reference_alleles'].append(allele)

        if allele_group not in allele_groups:
            allele_groups[allele_group] = []
        if allele not in allele_groups[allele_group]:
            allele_groups[allele_group].append(allele)

    output_file = f"output/processed_data/reference_alleles/{locus.lower()}.json"
    with open(output_file, 'w') as reference_alleles_file:
        alleles = json.dump(reference_alleles, reference_alleles_file)

    output_file = f"output/processed_data/allele_groups/{locus.lower()}.json"
    with open(output_file, 'w') as allele_groups_file:
        alleles = json.dump(allele_groups, allele_groups_file)


    print (reference_alleles)
    print ('')
    print (allele_groups)
