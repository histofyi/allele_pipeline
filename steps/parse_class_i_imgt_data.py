from typing import Dict
import json
from pathlib import Path

from common.allele import parse_hla_description, fasta_reader, process_sequence, allele_name_modifiers
from common.helpers import slugify

from rich import print

filename = f'tmp/ipd_imgt_hla_prot.fasta'


entries = []
truncated = []


pocket_residues = [7,9,24,45,59,62,63,66,67,69,70,73,74,76,77,80,81,84,95,97,99,114,116,118,143,147,150,152,156,158,159,163,167,171]


def generate_lists_for_locus(locus:str):

    protein_alleles = {}
    cytoplasmic_sequences = {}
    gdomain_sequences = {}
    pocket_pseudosequences = {}

    # counter for number of sequences in IMGT dataset
    imgt_count = 0
    for entry in fasta_reader(filename):
        allele_info = parse_hla_description(str(entry.description))

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


def parse_sequence_dict(sequences:Dict) -> Dict:
    for sequence in sequences:
        # set min to be far higher initially than the numerical representation of the allele number
        min = 10000000000
        for this_allele in sequences[sequence]['alleles']:
            allele_number = this_allele['gene_allele_name'].split('*')[1]
            if allele_number[-1] in allele_name_modifiers:
                allele_number = allele_number[:-1]
            while allele_number.count(':') < 3:
                    allele_number += ':01'
            allele_number_numeric = int(allele_number.replace(':',''))
            if allele_number_numeric < min:
                min = allele_number_numeric
                sequences[sequence]['canonical_allele'] = this_allele
    return sequences


def construct_class_i_hla_allele_lists(locus:str):

    imgt_count, protein_alleles, cytoplasmic_sequences, gdomain_sequences, pocket_pseudosequences = generate_lists_for_locus(locus)

    # now we'll iterate through the alleles to find the canonical allele (the one with the lowest number)
    for allele in protein_alleles:
        allele_count = len(protein_alleles[allele]['alleles'])
        canonical_sequence = None
        canonical_allele = None
        if allele_count > 1:
            min = 10000000000
            for this_allele in protein_alleles[allele]['alleles']:
                # first of all we'll split the gene name into what will become a string representation of the number
                allele_number = this_allele['gene_allele_name'].split('*')[1]
                # if there is a modifier on the end of the gene name, we'll remove it so we can cast to an int
                if allele_number[-1] in allele_name_modifiers:
                    allele_number = allele_number[:-1]
                # full length gene allele names have four components (three : ) so we'll pad any that are shorter
                while allele_number.count(':') < 3:
                    allele_number += ':01'
                # then cast to an int
                allele_number_numeric = int(allele_number.replace(':',''))
                # find the lowest numberic allele number
                if allele_number_numeric < min:
                    min = allele_number_numeric
                    canonical_allele = this_allele

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
            canonical_allele = protein_alleles[allele]['alleles'][0]
            canonical_sequence = protein_alleles[allele]['sequences'][0]
        protein_alleles[allele]['canonical_allele'] = canonical_allele
        protein_alleles[allele]['canonical_sequence'] = canonical_sequence



    for sequence_type in ['cytoplasmic_sequences', 'gdomain_sequences', 'pocket_pseudosequences']:
        directory_path = f"output/{sequence_type}"
        Path(directory_path).mkdir(parents=True, exist_ok=True)
        filename = f"{directory_path}/hla_{locus.lower()}.json"
        sequence_list = parse_sequence_dict(eval(sequence_type))
        with open(filename, "w") as json_file:
            json.dump(sequence_list, json_file, sort_keys=True, indent=4)

    directory_path = f"output/protein_alleles"
    Path(directory_path).mkdir(parents=True, exist_ok=True)
    filename = f"{directory_path}/hla_{locus.lower()}.json"
    with open(filename, "w") as json_file:
        json.dump(protein_alleles, json_file, sort_keys=True, indent=4)
    
    print ("")
    print (f"Total number of HLA-{locus} sequences: {imgt_count}")
    print (f"Total number of alleles: {len(protein_alleles)}")
    print (f"Total number of unique cytosolic sequences: {len(cytoplasmic_sequences)}")
    print (f"Total number of unique g-domain sequences: {len(gdomain_sequences)}")
    print (f"Total number of unique pocket pseudosequences: {len(pocket_pseudosequences)}")
    print (f"")

