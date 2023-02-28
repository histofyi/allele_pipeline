from typing import Dict, List
from Bio.SeqIO.FastaIO import FastaIterator

import hashlib


locus_ignore_list = ['MICA','MICB','TAP1','TAP2']

class_i_locus_list = ['B']

class_i_starts = [
    "GDTRPRY",
    "GFHSLRY",
    "GPHSLRY",
    "GSHSLWY",
    "GPHSMRY",
    "GSHSMRY",
    "GSHSLRY",
    "GSHSLKY",
    "VSHSMRY",
    "SSHSMRY",
    "RSHSMRY",
    "SSYSMRY",
    "GSHFLRY",
    "RSHSLQY",
    "AFHSLRY",
    "GSHSMKY",
    "GSHSLSY",
    "GYHSLNY",
    "SSHSLRY",
    "SHSFSRF",
    "GSHSLTY",
    "SSHSPRY",
    "GSHFMRY",
    "ASSHSLY",
    "GSHPMRY",
    "GSHSMSY",
    "DSHSMRY",
    "GSHSLGY",
    "CSHSLRY",
    "RSHSLRY",
    "GSRSMRY",
    "GLHSLRY",
    "RPHSLRY",
    "GPHSLSY",
    "GSHSMTY",
    "GSHSVRY",
    "GSHSLRF",
    "GFHSMRY",
    "ELHTLRY",
    "ELHSLRY",
    "EHHTLRY",
    "GWHSMRY",
    "GWHSLRY",
    "GSHSFSY",
    "GAHSMRY",
    "GSHCMRY",
    "GSHSLKH",
    "VSHSLRY",
    "DSHSMWY",
    "GSHSMRH",
    "CSHSMRY",
    "GSHWMRY",
    "GSHSLLY",
    "GTHSLRY",
    "GSHFMQY",
    "GSHSMWY",
    "GPRSLSY",
    "GLYSMRY",
    "GSQSMRY",
    "GSHSLLY",
    "GSHLFGL",
    "GLHSLNY",
    "GTHSLRY",
    "GSHFMQY",
    "GSHSMWY",
    "GSNSMRY",
    "ASHSMRY",
    "GSNSMRY",
    "GLPLLRY",
    "GSHSMRC",
    "CSHSMRD",
    "CSHSMKY",
    "CSRSMRY",
    "WSHSMRY",
    "CSHSMGY",
    "CFHSMRY",
    "CSHSVRY",
    "FSHSMRY",
    "YSHSMRY",
    "CSHSRRY",
    "CSHSTRY",
    "CSHCMRY",
    "CAHSMRY",
    "CSHSMRF",
    "CSHSMRC",
    "CSQSMRY",
    "CSHFMRY",
    "GSHSTRY",
    "CSPSMRY",
    "GSHSMGY",
    "GSHSISY",
    "SSHSMKY",
    "GSHSMGY",
    "GSHSKRY",
    "CSHSIRY",
    "GSPSMRY",
    "GCHSMRY",
    "GSLSMRY",
    "GSHSMRN",
    "GSLSMRY",
    "CSHSTKY",
    "GSDSMRY",
    "GTHSMRY",
    "GSYSMRY",
    "GSHTMRY",
    "GSHSRRY",
    "GSHSIRY",
    "GYHSMRY",
    "GSHSIRY",
    "GSHTMRY",
    "GSHSRRY",
    "GSHSIRY",
    "GSHAMRY",
    "GSHYMRY",
    "GSHAMRY",
    "GSHYMRY",
]

allele_name_modifiers = ['N','L','S','C','A','Q']

def fasta_reader(filename:str) -> List:
    with open(filename) as handle:
        for record in FastaIterator(handle):
            yield record


def parse_hla_description(description:str) -> Dict:
    description_elements = description.split(' ')
    # the locus is the part before the * character in the second element of the description
    locus = description_elements[1].split('*')[0]
    if locus not in locus_ignore_list:
        gene_allele_name = f'HLA-{description_elements[1]}'
    else:
        gene_allele_name = description_elements[1]
    protein_allele_name = ':'.join(gene_allele_name.split(":")[0:2])
    id = description_elements[0].split(':')[1]
    allele_info = {
        'protein_allele_name':protein_allele_name,
        'gene_allele_name':gene_allele_name,
        'locus': locus,
        'id':id
    }
    return allele_info


def process_sequence(sequence:str, locus:str) -> Dict:
    missing_start = None
    cytoplasmic_sequence = None
    appropriate_length = None
    start_match = False
    has_class_i_start = False
    if len(sequence) > 270:
        appropriate_length = True
        for start in class_i_starts:
            if start in sequence:
                has_class_i_start = True
                start_string = start
                missing_start = False
                start_match = True
                break
            elif start[1:] in sequence and not start_match:
                has_class_i_start = True
                start_string = start[1:]
                missing_start = True
    else:
        appropriate_length = False
    if has_class_i_start:
        cytoplasmic_sequence = start_string + sequence.split(start_string)[1]
        if missing_start:
            cytoplasmic_sequence = cytoplasmic_sequence[:274]
        else:
            cytoplasmic_sequence = cytoplasmic_sequence[:275]
        return {
            'cytoplasmic_sequence':cytoplasmic_sequence,
            'class_i_start':has_class_i_start,
            'missing_start':missing_start,
            'appropriate_length':appropriate_length
        }
    return {
            'cytoplasmic_sequence':cytoplasmic_sequence,
            'class_i_start':has_class_i_start,
            'missing_start':missing_start,
            'appropriate_length':appropriate_length
    }
