from common.allele import parse_hla_description, fasta_reader, process_sequence, class_i_locus_list, allele_name_modifiers
from common.helpers import slugify


filename = f'tmp/ipd_imgt_hla_prot.fasta'



entries = []
truncated = []

protein_alleles = {}

sequences = {}

# counter for number of sequences in IMGT dataset
imgt_count = 0
for entry in fasta_reader(filename):
    allele_info = parse_hla_description(str(entry.description))

    if allele_info['locus'] in class_i_locus_list:
        protein_allele_name = None
        allele_slug = None
        # increment the IMGT counter, only incrementing for class I sequences
        imgt_count += 1 
        sequence_data = process_sequence(str(entry.seq), allele_info['locus'])
        # we only want alleles with sequences long enough to include cytoplasmic domains
        if sequence_data['appropriate_length'] and sequence_data['cytoplasmic_sequence']: 
            if len(sequence_data['cytoplasmic_sequence']) > 270:
                protein_allele_name = allele_info['protein_allele_name']
                
                if protein_allele_name[-1] in allele_name_modifiers:
                    protein_allele_name = protein_allele_name[:-1]
                
                allele_slug = slugify(protein_allele_name)

                if allele_slug not in protein_alleles:
                    protein_alleles[allele_slug] = {
                        'sequences':[],
                        'alleles':[],
                        'canonical_allele':'',
                        'canonical_sequence':''
                    }
                protein_alleles[allele_slug]['alleles'].append(allele_info)
                if sequence_data['cytoplasmic_sequence'] not in  protein_alleles[allele_slug]['sequences']:
                     protein_alleles[allele_slug]['sequences'].append(sequence_data['cytoplasmic_sequence'])
                if sequence_data['cytoplasmic_sequence'] not in sequences:
                    sequences[sequence_data['cytoplasmic_sequence']] = {
                        'alleles':[],
                        'canonical_allele':{}
                    }
                sequences[sequence_data['cytoplasmic_sequence']]['alleles'].append(allele_info)




# now we'll iterate through the alleles to find the canonical allele (the one with the lowest number)
for allele in protein_alleles:
    #print (allele)
    allele_count = len(protein_alleles[allele]['alleles'])
    if allele_count > 1:
        canonical_allele = None
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
        canonical_allele = protein_alleles[allele]['alleles'][0]
        canonical_sequence = protein_alleles[allele]['sequences'][0]


for sequence in sequences:
    canonical_allele = None
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
    
    print (sequence)
    print (sequences[sequence]['canonical_allele'])
    print ('')

print (f"Total number of sequences: {imgt_count}")
print (f"Total number of alleles: {len(protein_alleles)}")
print (f"Total number of unique sequences: {len(sequences)}")


