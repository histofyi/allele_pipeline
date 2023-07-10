from typing import List, Dict, Union, Optional

netmhcpan_pseudosequence_file = f'tmp/MHC_pseudo.dat'


def parse_netmhcpan_allele_name(allele_name_string:str) -> Union[None, Dict]:
    allele_details = None
    if '-' in allele_name_string and (':' in allele_name_string or '*' in allele_name_string):
        allele_details = {}
        allele_details['species_stem'] = allele_name_string.split('-')[0]
        allele_info = allele_name_string.split('-')[1]
        if '*' in allele_info:
            allele_details['locus'] = allele_info.split('*')[0]
        else:
            if allele_details['species_stem'] == 'HLA' and allele_info[0] in ['A','B','C','E','F','G']:
                allele_details['locus'] = allele_info[0]
                allele_details['allele_group'] = f"{allele_details['species_stem']}-{allele_details['locus']}*{allele_info[1:].split(':')[0]}"
                allele_details['protein_allele_name'] = f"{allele_details['species_stem']}-{allele_details['locus']}*{allele_info[1:]}"
        allele_details['original_string'] = allele_name_string 
        allele_details['source'] = 'netmhcpan'
    return allele_details
    


with open(netmhcpan_pseudosequence_file, 'r') as netmhcpan:
    netmhcpan_data = netmhcpan.read()

pseudosequences = {}

i = 0
for row in netmhcpan_data.splitlines():
    components = row.split(' ')
    if len(components) > 1:
        allele = parse_netmhcpan_allele_name(components[0])
        pseudosequence = components[1]
        if allele is not None:
            if allele['species_stem'] == 'HLA':
                i += 1
                if pseudosequence not in pseudosequences:
                    pseudosequences[pseudosequence] = {'alleles':[]}
                else:
                    print ('already there')
                pseudosequences[pseudosequence]['alleles'].append(allele)   
                print (pseudosequence)
                print (allele)


print (i)
print (len(pseudosequences))