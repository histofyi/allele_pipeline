from typing import Dict, List, Tuple

import csv
import json


def get_allele_group(text:str) -> str:
    elements = text.split('_')
    return f"{elements[0]}_{elements[1]}_{elements[2]}".upper()



def deslugify_allele_group(text:str) -> str:
    elements = text.split('_')
    return f"{elements[0]}-{elements[1]}*{elements[2]}".upper()


def deslugify_allele(text:str) -> str:
    elements = text.split('_')
    return f"{elements[0]}-{elements[1]}*{elements[2]}:{elements[3]}".upper()


def deslugify_locus(test:str) -> str:
    elements = test.split('_')
    return f"{elements[0]}-{elements[1]}".upper()



def create_tabular_representations(config:Dict, **kwargs) -> None:
    """
    This function takes allele group lists for a locus and makes a tabular representation of them.

    Args:
        locus (str): the locus to be parsed
        verbose (bool): whether specific information is output to the terminal, for large sequence sets this can be overwhelming and significantly slow down the function
    Returns:

    """
    locus = kwargs['locus']
    species_stem = kwargs['species_stem']
    
    locus_slug = f"{species_stem}_{locus.lower()}"


    input_filename = f"output/processed_data/protein_alleles/{locus_slug}.json"
    output_filename = f"output/tabular_data/alleles/{locus_slug}.csv"

    with open(input_filename, 'r') as protein_alleles_filehandle:
        protein_alleles = json.load(protein_alleles_filehandle)
    
    pocket_positions = config['CONSTANTS']['IMGT_POCKET_RESIDUES']


    labels = ['allele_slug', 'allele', 'allele_id', 'allele_url', 'allele_group_slug', 'allele_group', 'locus_slug', 'locus', 'species', 'netmhcpan_pseudosequence']

    for position in pocket_positions:
        labels.append(f"alpha_{position}")

    table = []
    table.append(labels)

    species_slug = 'homo_sapiens'
    for allele_slug in protein_alleles:
        allele_id = protein_alleles[allele_slug]['alleles'][0]['id']
        row = [
            allele_slug,
            deslugify_allele(allele_slug),
            f"imgt/hla:{allele_id}",
            f"https://www.ebi.ac.uk/ipd/imgt/hla/alleles/allele/?accession={allele_id}",
            get_allele_group(allele_slug),
            deslugify_allele_group(allele_slug),
            locus,
            deslugify_locus(locus_slug),
            species_slug,
            protein_alleles[allele_slug]['pocket_pseudosequence'],
        ]
        for position in protein_alleles[allele_slug]['pocket_pseudosequence']:
            row.append(position)
        table.append(row)
        

    with open(output_filename, 'w', newline='\n') as f:
        writer = csv.writer(f)
        writer.writerows(table)



    pass
