from typing import Dict
import json


from common.allele import parse_hla_description, parse_h2_description, fasta_reader, process_sequence, find_canonical_allele, allele_name_modifiers
from common.helpers import slugify

from rich import print



entries = []
truncated = []

species_functions = {
    'hla': parse_hla_description,
    'h2': parse_h2_description
}

pocket_residues = [7,9,24,45,59,62,63,66,67,69,70,73,74,76,77,80,81,84,95,97,99,114,116,118,143,147,150,152,156,158,159,163,167,171]

non_standard_nomenclature_species = ['h2']

def generate_lists_for_locus(locus:str, species_slug, sequence_set:str, verbose:bool):
    """
    This function takes a dataset and generate an allele list and associated sequence lists for a specific locus.

    Args:
        locus (str): the locus of interest e.g. A
        species_slug (str): the slug for the species, this is used to switch between allele numbering functions for the mouse in particular e.g. h2, hla
        sequence_set (str): the name of the sequence set, this is used to determine the file name e.g. IPD_IMGT_HLA_PROT which results in the filename tmp/ipd_imgt_hla_prot.fasta
        verbose (bool): whether specific information is output to the terminal, for large sequence sets this can be overwhelming and significantly slow down the function
    Returns:

    """
    filename = f"tmp/{sequence_set.lower()}.fasta"

    protein_alleles = {}
    cytoplasmic_sequences = {}
    gdomain_sequences = {}
    pocket_pseudosequences = {}

    # counter for number of sequences in IMGT dataset
    imgt_count = 0
    for entry in fasta_reader(filename):
        if verbose:
            print (entry)
            print (species_slug)
        
        allele_info = species_functions[species_slug](str(entry.description))

        if verbose:
            print (allele_info)

        if allele_info['locus'] in locus:
            protein_allele_name = None
            allele_slug = None
            # increment the IMGT counter, only incrementing for class I sequences
            imgt_count += 1 
            sequence_data = process_sequence(str(entry.seq), pocket_residues)
            # we only want alleles with sequences long enough to include cytoplasmic domains
            if sequence_data['appropriate_length'] and sequence_data['cytoplasmic_sequence']: 
                if len(sequence_data['cytoplasmic_sequence']) > 270:
                    protein_allele_name = allele_info['protein_allele_name']
                    
                    if species_slug not in non_standard_nomenclature_species:
                        # we need to check if the allele name has a modifier, such as N, Q etc
                        if protein_allele_name[-1] in allele_name_modifiers:
                            protein_allele_name = protein_allele_name[:-1]
                    
                    # slugify the cleaned allele name
                    allele_slug = slugify(protein_allele_name)

                    # and check if it's in the protein allele dict
                    if allele_slug not in protein_alleles:
                        protein_alleles[allele_slug] = {
                            'sequences':[],
                            'alleles':[],
                            'canonical_allele':'',
                            'canonical_sequence':'',
                            'gdomain_sequence':sequence_data['gdomain_sequence'],
                            'pocket_pseudosequence':sequence_data['pocket_pseudosequence']
                        }
                    # and append the specific allele information    
                    protein_alleles[allele_slug]['alleles'].append(allele_info)

                    if sequence_data['cytoplasmic_sequence'] not in  protein_alleles[allele_slug]['sequences']:
                        protein_alleles[allele_slug]['sequences'].append(sequence_data['cytoplasmic_sequence'])
                    
                    # now add the unique sequences to the different sequence dictionaries
                    # first up the cytoplasmic sequence dict 
                    if sequence_data['cytoplasmic_sequence'] not in cytoplasmic_sequences:
                        cytoplasmic_sequences[sequence_data['cytoplasmic_sequence']] = {
                            'alleles':[],
                            'canonical_allele':{}
                        }
                    cytoplasmic_sequences[sequence_data['cytoplasmic_sequence']]['alleles'].append(allele_info)
                    
                    # next the gdomain dict
                    if sequence_data['gdomain_sequence'] not in gdomain_sequences:
                        gdomain_sequences[sequence_data['gdomain_sequence']] = {
                            'alleles':[],
                            'canonical_allele':{}
                        }
                    gdomain_sequences[sequence_data['gdomain_sequence']]['alleles'].append(allele_info)
                    
                    # finally the pocket pseudosequence dict (like NetMHCPan)
                    if sequence_data['pocket_pseudosequence'] not in pocket_pseudosequences:
                        pocket_pseudosequences[sequence_data['pocket_pseudosequence']] = {
                            'alleles':[],
                            'canonical_allele':{}
                        }
                    pocket_pseudosequences[sequence_data['pocket_pseudosequence']]['alleles'].append(allele_info)
    
    return imgt_count, protein_alleles, cytoplasmic_sequences, gdomain_sequences, pocket_pseudosequences


def parse_sequence_dict(sequences:Dict, species_slug:str) -> Dict:
    """
    This function takes a dictionary of processed sequences (for example gdomain sequences) and finds the canonical allele, assuming that is the one with the lowest allele number
    
    Args:

    """
    if species_slug not in non_standard_nomenclature_species:
        for sequence in sequences:
            sequences[sequence] = find_canonical_allele(sequences[sequence])
    else:
        # TODO think about whether this rather blunt approach is the correct one for species such as the mouse with non standard allele nomenclature
        # TODO think about wehther this should be handled in find_canonical_allele instead
        for sequence in sequences:
            sequences[sequence]['canonical_allele'] = sequences[sequence]['alleles'][0]
    return sequences


def construct_class_i_locus_allele_lists(locus:str, species_slug:str, sequence_set:str, verbose=False):

    imgt_count, protein_alleles, cytoplasmic_sequences, gdomain_sequences, pocket_pseudosequences = generate_lists_for_locus(locus, species_slug, sequence_set, verbose)

    # now we'll iterate through the alleles to find the canonical allele (the one with the lowest number)
    for allele in protein_alleles:
        allele_count = len(protein_alleles[allele]['alleles'])
        canonical_sequence = None
        canonical_allele = None
        if allele_count > 1 and species_slug not in non_standard_nomenclature_species:
            protein_alleles[allele] = find_canonical_allele(protein_alleles[allele])

            # there may be many different length variants of the sequence, we just want the longest one to be the canonical one for matching
            if len(protein_alleles[allele]['sequences']) > 1:
                max_sequence_length = 0
                canonical_sequence = ''
                for sequence in protein_alleles[allele]['sequences']:
                    if len(sequence) > max_sequence_length:
                        max_sequence_length = len(sequence)
                        canonical_sequence = sequence
            else:
                canonical_sequence = protein_alleles[allele]['sequences'][0]
        else:
            protein_alleles[allele]['canonical_allele'] = protein_alleles[allele]['alleles'][0]
            canonical_sequence = protein_alleles[allele]['sequences'][0]
        protein_alleles[allele]['canonical_sequence'] = canonical_sequence



    for sequence_type in ['cytoplasmic_sequences', 'gdomain_sequences', 'pocket_pseudosequences']:
        # generate a directory path for each type of sequence
        directory_path = f"output/{sequence_type}"
        # create the directory path if it doesn't exist
        # generate a filename for the specific file for the locus
        filename = f"{directory_path}/{species_slug}_{locus.lower()}.json"
        # parse the sequence dictionary to assign the canonical allele
        sequence_list = parse_sequence_dict(eval(sequence_type), species_slug)
        # write the sequence dictionary
        with open(filename, "w") as json_file:
            json.dump(sequence_list, json_file, sort_keys=True, indent=4)

    # now generate a directory path 
    directory_path = f"output/protein_alleles"
    filename = f"{directory_path}/{species_slug}_{locus.lower()}.json"
    with open(filename, "w") as json_file:
        json.dump(protein_alleles, json_file, sort_keys=True, indent=4)
    
    print ("")
    print (f"Total number of {species_slug.upper()}-{locus} sequences: {imgt_count}")
    print (f"Total number of alleles: {len(protein_alleles)}")
    print (f"Total number of unique cytosolic sequences: {len(cytoplasmic_sequences)}")
    print (f"Total number of unique g-domain sequences: {len(gdomain_sequences)}")
    print (f"Total number of unique pocket pseudosequences: {len(pocket_pseudosequences)}")
    print (f"")

