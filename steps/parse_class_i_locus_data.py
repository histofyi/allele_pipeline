from typing import Dict, Union
import json

from common.allele import parse_hla_description, parse_h2_description, fasta_reader, process_sequence, find_canonical_allele, allele_name_modifiers
from common.helpers import slugify

from rich import print


species_functions = {
    'hla': parse_hla_description,
    'h2': parse_h2_description
}

non_standard_nomenclature_species = ['h2']

def generate_lists_for_locus(locus:str, species_slug, sequence_set:str, config:Dict, verbose:bool) -> Union[int, Dict, Dict, Dict, Dict]:
    """
    This function takes a dataset and generate an allele list and associated sequence lists for a specific locus.

    Args:
        locus (str): the locus of interest e.g. A
        species_slug (str): the slug for the species, this is used to switch between allele numbering functions for the mouse in particular e.g. h2, hla
        sequence_set (str): the name of the sequence set, this is used to determine the file name e.g. IPD_IMGT_HLA_PROT which results in the filename tmp/ipd_imgt_hla_prot.fasta
        verbose (bool): whether specific information is output to the terminal, for large sequence sets this can be overwhelming and significantly slow down the function
    Returns:
        int: the number of class I sequences within a specific locus in the dataset
        Dict: the dictionary of protein alleles 
        Dict: the dictionary of cytoplasmic sequences
        Dict: the dictionary of g-domain sequences
        Dict: the dictionary of pocket pseudosequences (same as NetMHCPan pseudosequences)
    """
    filename = f"tmp/{sequence_set.lower()}.fasta"

    protein_alleles = {}
    cytoplasmic_sequences = {}
    gdomain_sequences = {}
    pocket_pseudosequences = {}

    # counter for number of sequences in the dataset
    total_sequences = 0
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
            # increment the total sequence counter, only incrementing for class I sequences
            total_sequences += 1 
            sequence_data = process_sequence(str(entry.seq), config['CONSTANTS']['IMGT_POCKET_RESIDUES'])
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
                    # TODO this looks optimisable

                    for sequence_type in config['CONSTANTS']['SEQUENCE_TYPES']:
                        
                        this_sequence_type = eval(sequence_type)
                        # sequence types in the config are plural, in the protein_allele dictionary they're singular
                        sequence_type = sequence_type[:-1]

                        # if the sequence is not in the specific sequence type dictionary then we need to create an entry
                        if sequence_data[sequence_type] not in this_sequence_type:
                            this_sequence_type[sequence_data[sequence_type]] = {
                                'alleles':[],
                                'canonical_allele':{}
                            }
                        # and then append the allele info to the alleles list for the sequence
                        this_sequence_type[sequence_data[sequence_type]]['alleles'].append(allele_info)

    return total_sequences, protein_alleles, cytoplasmic_sequences, gdomain_sequences, pocket_pseudosequences


def parse_sequence_dict(sequences:Dict, species_slug:str) -> Dict:
    """
    This function takes a dictionary of processed sequences (for example gdomain sequences) and finds and appends the canonical allele, assuming that is the one with the lowest allele number
    
    Args:
        sequences (Dict): a dictionary of sequences
        species_slug (str): a slugified version of the species prefix e.g. hla, h2, gaga

    Returns:
        Dict: the marked up version of the sequences dictionary
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


def construct_class_i_locus_allele_lists(config:Dict, **kwargs) -> Dict:
    """
    This function creates the data for a specific locus and saves it in a set of files in the output directory

    Args:
        locus (str): the locus to be processed in a dataset e.g. A for HLA-A
        species_slug (str): the slugified version of the species stem for the sequences e.g. hla for HLA
        sequence_set (str): the dataset to be processed e.g. IPD_IMGT_HLA_PROT for the HLA protein sequence dataset from IPD/IMGT
        config (Dict): the configuration dictionary
        verbose (bool): a boolean as to whether this step should output to the terminal

    Returns:
        Dict: the action dictionary for this step which will be stored in the pipeline log
    """
    
    locus = kwargs['locus']
    species_slug = kwargs['species_slug']
    sequence_set = kwargs['sequence_set']
    output_path = kwargs['output_path']
    if 'verbose' in kwargs:
        verbose = kwargs['verbose']
    else:
        verbose = False

    total_sequences, protein_alleles, cytoplasmic_sequences, gdomain_sequences, pocket_pseudosequences = generate_lists_for_locus(locus, species_slug, sequence_set, config, verbose)

    # now we'll iterate through the alleles to find the canonical allele (the one with the lowest number)
    for allele in protein_alleles:
        allele_count = len(protein_alleles[allele]['alleles'])
        canonical_sequence = None
        canonical_allele = None

        # if there's more than one sequence and the species is not one with non-standard nomenclature, look for the canonical allele and sequence
        if allele_count > 1 and species_slug not in non_standard_nomenclature_species:
            protein_alleles[allele] = find_canonical_allele(protein_alleles[allele])

            # there may be many different length variants of the sequence, we just want the longest one to be the canonical one for matching
            if len(protein_alleles[allele]['sequences']) > 1:
                max_sequence_length = 0
                canonical_sequence = ''
                # now iterate through the sequences
                for sequence in protein_alleles[allele]['sequences']:
                    # find the longest sequence, and make it the canonical one
                    if len(sequence) > max_sequence_length:
                        max_sequence_length = len(sequence)
                        canonical_sequence = sequence
            else:
                # or if there's no sequence variation, just make the first sequence in the array the canonical one
                canonical_sequence = protein_alleles[allele]['sequences'][0]
        else:
            # if there is only one allele, the only allele is the canonical one and canonical sequence
            protein_alleles[allele]['canonical_allele'] = protein_alleles[allele]['alleles'][0]
            canonical_sequence = protein_alleles[allele]['sequences'][0]
        protein_alleles[allele]['canonical_sequence'] = canonical_sequence


    # next we'll iterate through each type of sequence and process them and save each one to a file
    for sequence_type in ['cytoplasmic_sequences', 'gdomain_sequences', 'pocket_pseudosequences']:
        
        # generate a directory path for each type of sequence
        directory_path = f"{output_path}/processed_data/{sequence_type}"

        # generate a filename for the specific file for the locus
        filename = f"{directory_path}/{species_slug}_{locus.lower()}.json"
        
        # parse the sequence dictionary to assign the canonical allele
        sequence_list = parse_sequence_dict(eval(sequence_type), species_slug)
        
        # write the sequence dictionary
        with open(filename, "w") as json_file:
            json.dump(sequence_list, json_file, sort_keys=True, indent=4)

    # now generate a dictionary file for alleles 
    directory_path = f"{output_path}/processed_data/protein_alleles"
    filename = f"{directory_path}/{species_slug}_{locus.lower()}.json"
    with open(filename, "w") as json_file:
        json.dump(protein_alleles, json_file, sort_keys=True, indent=4)

    # output some statistics to the terminal if verbose is True
    if verbose:
        print ("")
        print (f"Total number of {species_slug.upper()}-{locus} sequences: {total_sequences}")
        print (f"Total number of alleles: {len(protein_alleles)}")
        print (f"Total number of unique cytoplasmic sequences: {len(cytoplasmic_sequences)}")
        print (f"Total number of unique g-domain sequences: {len(gdomain_sequences)}")
        print (f"Total number of unique pocket pseudosequences: {len(pocket_pseudosequences)}")
        print (f"")
    
    action_log = {
        'locus': f"{species_slug.upper()}-{locus}",
        'sequences_processed': total_sequences,
        'alleles_found': len(protein_alleles),
        'unique_cytoplasmic_sequences': len(cytoplasmic_sequences),
        'unique_gdomain_sequences': len(gdomain_sequences),
        'unique_pocket_pseudosequences': len(pocket_pseudosequences)
    }

    return action_log
