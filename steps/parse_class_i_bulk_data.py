from typing import Dict, List, Union
import json

from common.allele import parse_mhc_description, fasta_reader, process_sequence, find_canonical_allele, allele_name_modifiers
from common.helpers import slugify

from rich import print


# TODO handle non-IMGT conforming sequences better at the moment this is set to None so that incorrect pocket pseuedosequences aren't generated
pocket_residues = None


def generate_lists(sequence_set:str, verbose:bool=True) -> Union[Dict, Dict, Dict, Dict, Dict, List]:
    """
    This function takes a dataset and generate an allele list and associated sequence lists for all Class I loci contained within it.

    Args:
        sequence_set (str): the name of the sequence set, this is used to determine the file name e.g. IPD_MHC_PROT which results in the filename tmp/ipd_mhc_prot.fasta
        verbose (bool): whether specific information is output to the terminal, for large sequence sets this can be overwhelming and significantly slow down the function
    Returns:
        Dict: the dictionary of protein alleles 
        Dict: the dictionary of cytoplasmic sequences
        Dict: the dictionary of g-domain sequences
        Dict: the dictionary of pocket pseudosequences (same as NetMHCPan pseudosequences), in this case empty until we can generate them reliably for non-IMGT numbering compliant sequences
        Dict: the dictionary of basic statistics
        List: an array of unmatched alleles
    """
    # this variable stores a simple counter of all sequences in the dataset
    all_sequences_count = 0
    
    filename = f"tmp/{sequence_set.lower()}.fasta"

    unmatched = []

    protein_alleles = {}
    cytoplasmic_sequences = {}
    gdomain_sequences = {}
    pocket_pseudosequences = None

    mhc_class_i_count = 0
    for entry in fasta_reader(filename):
        allele_info = parse_mhc_description(str(entry.description))
        
        # first of all exclude any sequences whose locus is that of classical Class II, this won't catch all species, but will do for primates
        if not allele_info['locus'].split('-')[1][0:3] in ['DPA','DPB','DQA','DQB','DRA','DRB']:
            
            allele_slug = None
            locus_slug = None

            sequence_data = process_sequence(str(entry.seq), pocket_residues)
            if sequence_data['class_i_start'] and sequence_data['appropriate_length']:
                mhc_class_i_count += 1

                # slugify the allele name
                allele_slug = slugify(allele_info['protein_allele_name'])

                # slugify the locus
                locus_slug = slugify(allele_info['locus'])
                
                # slugify the locus
                species_slug = slugify(allele_info['locus'].split('-')[0])

                if locus_slug not in protein_alleles:
                    protein_alleles[locus_slug] = {}

                # and check if it's in the protein allele dict
                if allele_slug not in protein_alleles[locus_slug]:
                    protein_alleles[locus_slug][allele_slug] = {
                        'sequences':[],
                        'alleles':[],
                        'canonical_allele':'',
                        'canonical_sequence':'',
                        'gdomain_sequence':sequence_data['gdomain_sequence'],
                        'pocket_pseudosequence':sequence_data['pocket_pseudosequence']
                    }
                # and append the specific allele information    
                protein_alleles[locus_slug][allele_slug]['alleles'].append(allele_info)

                if sequence_data['cytoplasmic_sequence'] not in  protein_alleles[locus_slug][allele_slug]['sequences']:
                    protein_alleles[locus_slug][allele_slug]['sequences'].append(sequence_data['cytoplasmic_sequence'])
                # now add the unique sequences to the different sequence dictionaries
                # first up the cytoplasmic sequence dict 
                if not species_slug in cytoplasmic_sequences:
                    cytoplasmic_sequences[species_slug] = {}
                if not locus_slug in cytoplasmic_sequences[species_slug]:
                    cytoplasmic_sequences[species_slug][locus_slug] = {}

                if sequence_data['cytoplasmic_sequence'] not in cytoplasmic_sequences[species_slug][locus_slug]:
                    cytoplasmic_sequences[species_slug][locus_slug][sequence_data['cytoplasmic_sequence']] = {
                        'alleles':[],
                        'canonical_allele':{}
                    }
                cytoplasmic_sequences[species_slug][locus_slug][sequence_data['cytoplasmic_sequence']]['alleles'].append(allele_info)
                
                # next the gdomain dict
                if not species_slug in gdomain_sequences:
                    gdomain_sequences[species_slug] = {}
                if not locus_slug in gdomain_sequences[species_slug]:
                    gdomain_sequences[species_slug][locus_slug] = {}

                if sequence_data['gdomain_sequence'] not in gdomain_sequences[species_slug][locus_slug]:
                    gdomain_sequences[species_slug][locus_slug][sequence_data['gdomain_sequence']] = {
                        'alleles':[],
                        'canonical_allele':{}
                    }
                gdomain_sequences[species_slug][locus_slug][sequence_data['gdomain_sequence']]['alleles'].append(allele_info)
            else:
                unmatched.append(allele_info['protein_allele_name'])
                
        all_sequences_count += 1
    stats = {
        'all_sequences_count':all_sequences_count,
        'all_class_i_allele_sequences_count':mhc_class_i_count,
        'unmatched_count': len(unmatched)
    }
    return protein_alleles, cytoplasmic_sequences, gdomain_sequences, pocket_pseudosequences, stats, unmatched


def parse_sequence_dict(sequences:Dict) -> Dict:
    for sequence in sequences:
        sequences[sequence] = find_canonical_allele(sequences[sequence])
    return sequences


def construct_class_i_bulk_allele_lists(config:Dict, **kwargs):
    """
    This function creates the data for a specific locus and saves it in a set of files in the output directory

    Args:
        config (Dict): the configuration dictionary
        sequence_set (str): the dataset to be processed e.g. IPD_IMGT_HLA_PROT for the HLA protein sequence dataset from IPD/IMGT
        verbose (bool): a boolean as to whether this step should output to the terminal

    Returns:
        Dict: the action dictionary for this step which will be stored in the pipeline log
    """
    sequence_set = kwargs['sequence_set']
    output_path = kwargs['output_path']
    if 'verbose' in kwargs:
        verbose = kwargs['verbose']
    else:
        verbose = False
    
    protein_alleles, cytoplasmic_sequences, gdomain_sequences, pocket_pseudosequences, stats, unmatched = generate_lists(sequence_set)

    all_allele_count = 0
    # now we'll iterate through the alleles to find the canonical allele (the one with the lowest number)
    for locus_slug in protein_alleles:
        for allele in protein_alleles[locus_slug]:
            all_allele_count +=1
            allele_count = len(protein_alleles[locus_slug][allele]['alleles'])
            canonical_sequence = None
            if allele_count > 1:
                protein_alleles[locus_slug][allele] = find_canonical_allele(protein_alleles[locus_slug][allele])
                
                # there may be many different length variants of the sequence, we just want the longest one to be the canonical one for matching
                if len(protein_alleles[locus_slug][allele]['sequences']) > 1:
                    max_sequence_length = 0
                    canonical_sequence = ''
                    for sequence in protein_alleles[locus_slug][allele]['sequences']:
                        if len(sequence) > max_sequence_length:
                            max_sequence_length = len(sequence)
                            canonical_sequence = sequence
                else:
                    canonical_sequence = protein_alleles[locus_slug][allele]['sequences'][0]
            else:
                protein_alleles[locus_slug][allele]['canonical_allele'] = protein_alleles[locus_slug][allele]['alleles'][0]
                canonical_sequence = protein_alleles[locus_slug][allele]['sequences'][0]
            protein_alleles[locus_slug][allele]['canonical_sequence'] = canonical_sequence


    for sequence_type in ['cytoplasmic_sequences', 'gdomain_sequences', 'pocket_pseudosequences']:
        directory_path = f"{output_path}/processed_data/{sequence_type}"
        sequences = eval(sequence_type)
        if sequences:
            for species in sequences:
                for locus in sequences[species]:
                    filename = f"{directory_path}/{locus.lower()}.json"
                    sequence_list = parse_sequence_dict(sequences[species][locus])
                    with open(filename, "w") as json_file:
                        json.dump(sequence_list, json_file, sort_keys=True, indent=4)
    
    directory_path = f"{output_path}/processed_data/protein_alleles"
    for locus_slug in protein_alleles:
        filename = f"{directory_path}/{locus_slug}.json"
        with open(filename, "w") as json_file:
            json.dump(protein_alleles[locus_slug], json_file, sort_keys=True, indent=4)

    action_log = {k:v for k,v in stats.items()}

    action_log['species_found'] = len(cytoplasmic_sequences)
    action_log['loci_found'] = len(protein_alleles)
    action_log['unique_class_i_alleles_found'] = all_allele_count
    
    return action_log


