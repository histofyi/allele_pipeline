from typing import Dict
import json

from common.allele import parse_mhc_description, fasta_reader, process_sequence, allele_name_modifiers
from common.helpers import slugify

from rich import print

filename = f'tmp/ipd_mhc_prot.fasta'


entries = []
truncated = []


pocket_residues = None



def generate_lists():
    i = 0

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
                unmatched.append(allele_info['locus'])
                
        i += 1

    print (f"MHC Class I sequences found: {mhc_class_i_count}")
    print (f"Number of unmatched sequences to analyse: {len(unmatched)}")   
    
    return i, protein_alleles, cytoplasmic_sequences, gdomain_sequences, pocket_pseudosequences


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

            if allele_number.replace(':','').isnumeric():
                # then cast to an int
                allele_number_numeric = int(allele_number.replace(':',''))
                # find the lowest numberic allele number
                if allele_number_numeric < min:
                    min = allele_number_numeric
                    sequences[sequence]['canonical_allele'] = this_allele
    return sequences


def construct_class_i_bulk_allele_lists(dataset_name:str):
    ipd_mhc_count, protein_alleles, cytoplasmic_sequences, gdomain_sequences, pocket_pseudosequences = generate_lists()

    # now we'll iterate through the alleles to find the canonical allele (the one with the lowest number)
    for locus_slug in protein_alleles:
        for allele in protein_alleles[locus_slug]:
            allele_count = len(protein_alleles[locus_slug][allele]['alleles'])
            canonical_sequence = None
            canonical_allele = None
            if allele_count > 1:
                min = 10000000000
                for this_allele in protein_alleles[locus_slug][allele]['alleles']:
                    # first of all we'll split the gene name into what will become a string representation of the number
                    allele_number = this_allele['gene_allele_name'].split('*')[1]
                    # if there is a modifier on the end of the gene name, we'll remove it so we can cast to an int
                    if allele_number[-1] in allele_name_modifiers:
                        allele_number = allele_number[:-1]
                    # full length gene allele names have four components (three : ) so we'll pad any that are shorter
                    while allele_number.count(':') < 3:
                        allele_number += ':01'

                    if allele_number.replace(':','').isnumeric():
                        # then cast to an int
                        allele_number_numeric = int(allele_number.replace(':',''))
                        # find the lowest numberic allele number
                        if allele_number_numeric < min:
                            min = allele_number_numeric
                            canonical_allele = this_allele
                
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
                canonical_allele = protein_alleles[locus_slug][allele]['alleles'][0]
                canonical_sequence = protein_alleles[locus_slug][allele]['sequences'][0]
            protein_alleles[locus_slug][allele]['canonical_allele'] = canonical_allele
            protein_alleles[locus_slug][allele]['canonical_sequence'] = canonical_sequence


    for sequence_type in ['cytoplasmic_sequences', 'gdomain_sequences', 'pocket_pseudosequences']:
        directory_path = f"output/{sequence_type}"
        sequences = eval(sequence_type)
        if sequences:
            for species in sequences:
                for locus in sequences[species]:
                    filename = f"{directory_path}/{locus.lower()}.json"
                    sequence_list = parse_sequence_dict(sequences[species][locus])
                    with open(filename, "w") as json_file:
                        json.dump(sequence_list, json_file, sort_keys=True, indent=4)
    
    directory_path = f"output/protein_alleles"
    for locus_slug in protein_alleles:
        filename = f"{directory_path}/{locus_slug}.json"
        with open(filename, "w") as json_file:
            json.dump(protein_alleles[locus_slug], json_file, sort_keys=True, indent=4)


    print (f"Total sequences in {dataset_name} dataset: {ipd_mhc_count}")
    
    print (f"Number of species found: {len(cytoplasmic_sequences)}")
    print (f"Number of loci found: {len(protein_alleles)}")
    


