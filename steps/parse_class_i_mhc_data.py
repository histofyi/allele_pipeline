from typing import Dict
import json
from pathlib import Path

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
                

                if locus_slug not in protein_alleles:
                    protein_alleles[locus_slug] = {}
                    print (locus_slug)

                # and check if it's in the protein allele dict
                if allele_slug not in protein_alleles[locus_slug]:
                    protein_alleles[locus_slug][allele_slug] = {
                        'sequences':[],
                        'alleles':[],
                        'canonical_allele':'',
                        'canonical_sequence':''
                    }
                # and append the specific allele information    
                protein_alleles[locus_slug][allele_slug]['alleles'].append(allele_info)

            else:
                unmatched.append(allele_info['locus'])
                
        i += 1

    print (i)
    print (mhc_class_i_count)
    print (len(unmatched))
    print (len(protein_alleles))

    for locus in protein_alleles:
        print (locus.upper())
        for allele in protein_alleles[locus]:
            print (allele)
            print (len(protein_alleles[locus][allele]['alleles']))
    pass


def construct_class_i_mhc_allele_lists():

    
    pass


generate_lists()