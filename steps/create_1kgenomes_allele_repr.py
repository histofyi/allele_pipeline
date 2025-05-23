from common.onekgenomes import read_onekgenomes_data, process_alleles_and_groups

import json




def create_1kgenomes_json():
    headers, rows = read_onekgenomes_data()

    raw_alleles, raw_allele_groups, sample_count = process_alleles_and_groups(headers, rows, 'class_i')

    output_folder = 'output/processed_data/1kgenomes'

    loci = ['hla_a', 'hla_b', 'hla_c']

    alleles = {}
    allele_groups = {}

    for allele_group in raw_allele_groups:
        for locus in loci:
            if locus in allele_group:
                if locus not in allele_groups:
                    allele_groups[locus] = {}
                allele_groups[locus][allele_group] = raw_allele_groups[allele_group]

    for allele in raw_alleles:
        for allele_group in raw_allele_groups:
            if allele_group in allele:
                if allele_group not in alleles:
                    alleles[allele_group] = {}
                alleles[allele_group][allele] = raw_alleles[allele]
                

    with open(f"{output_folder}/1k_alleles.json", 'w') as f:
        json.dump(alleles, f, indent=4)

    with open(f"{output_folder}/1k_allele_groups.json", 'w') as f:
        json.dump(allele_groups, f, indent=4)



if __name__ == '__main__':
    create_1kgenomes_json()