import csv
import json

from .helpers import slugify


class_i_fields = ['HLA-A 1', 'HLA-A 2', 'HLA-B 1', 'HLA-B 2', 'HLA-C 1', 'HLA-C 2']
class_ii_fields = ['HLA-DQB1 1', 'HLA-DQB1 2', 'HLA-DRB1 1', 'HLA-DRB1 2']

sample_id_field = 'Sample ID'

super_population_field = 'Region'
population_field = 'Population'



def read_onekgenomes_data():
    headers = []
    rows = []
    with open('data/1000_genomes_hla.txt', 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                headers = [header for header in line.strip().split('\t')]
            else:
                row = [value for value in line.strip().split('\t')]
                rows.append(row)
    return headers, rows


def cleanup_allele_numbers(value):
    if '/' in value:
        value = value.split('/')[0]
    if ' ' in value:
        value = value.split(' ')[0]
    if value.count('*') == 1:
        value = value.split('*')[0]
    if 'N' in value:
        value = None
    return value
    

def update_dictionary(dictionary, key, super_population, population, dictionary_type):
    slug = slugify(key)
    if slug not in dictionary:
        dictionary[slug] = {'populations': {}, 'super_populations': {}, dictionary_type: key, 'count': 0}
    dictionary[slug]['count'] += 1

    if population not in dictionary[slug]['populations']:
        dictionary[slug]['populations'][population] = 0
    dictionary[slug]['populations'][population] += 1
                
    if super_population not in dictionary[slug]['super_populations']:
        dictionary[slug]['super_populations'][super_population] = 0
    dictionary[slug]['super_populations'][super_population] += 1

    return dictionary


def update_simple_dictionary(dictionary, key):
    if key not in dictionary:
        dictionary[key] = 0
    dictionary[key] += 1
    return dictionary


def min_max_normalise(value, min_value, max_value):
    if value != 0 and max_value:
        if min_value == max_value:
            return None
        else:    
            return round((value - min_value) / (max_value - min_value), 3)
    else:
        return None


def percentage(value, total):
    if total != 0:
        return value * 100 / total
    else:
        return None


def post_process_dictionary(dictionary, populations, super_populations):
    raw_values = {}
    processed_dictionary = {}

    keys = sorted(list(dictionary.keys()))

    loci = {}

    for key in keys:
        for locus in ['hla_a', 'hla_b', 'hla_c']:
            if locus in key:
                if locus not in loci:
                    loci[locus] = []
                loci[locus].append(key)

    raw_values['ALL'] = {}

    collections = list(super_populations.keys())


    for locus in loci:
        keys = loci[locus]
        raw_values['ALL'][locus] = {}
        raw_values['ALL'][locus]['counts'] = [dictionary[key]['count'] for key in keys]
        raw_values['ALL'][locus]['min'] = min(raw_values['ALL'][locus]['counts'])
        raw_values['ALL'][locus]['max'] = max(raw_values['ALL'][locus]['counts'])
        raw_values['ALL'][locus]['total'] = sum(super_populations.values())

    
        for super_population in super_populations:
            if super_population not in raw_values:
                raw_values[super_population] = {}
            raw_values[super_population][locus] = {}
            raw_values[super_population][locus]['counts'] = [dictionary[key]['super_populations'][super_population] if super_population in dictionary[key]['super_populations'] else 0 for key in keys]
            raw_values[super_population][locus]['total'] = sum(raw_values[super_population][locus]['counts'])
            raw_values[super_population][locus]['percentages'] = [percentage(dictionary[key]['super_populations'][super_population], raw_values[super_population][locus]['total']) if super_population in dictionary[key]['super_populations'] else 0 for key in keys]
            raw_values[super_population][locus]['min'] = min(raw_values[super_population][locus]['percentages'])
            raw_values[super_population][locus]['max'] = max(raw_values[super_population][locus]['percentages'])
            


    print (collections)


    for locus in loci:
        keys = loci[locus]
        for i, key in enumerate(keys):
            if key not in processed_dictionary:
                processed_dictionary[key] = {}
            for collection in collections:
                if collection not in processed_dictionary[key]:
                    processed_dictionary[key][collection] = {}
                count = raw_values[collection][locus]['counts'][i]
                percentage_value = raw_values[collection][locus]['percentages'][i]
                total = raw_values[collection][locus]['total']
                min_value = raw_values[collection][locus]['min']
                max_value = raw_values[collection][locus]['max']
                min_max_normalised = min_max_normalise(percentage_value, min_value, max_value)
                processed_dictionary[key][collection] = {
                    'count': count,
                    'percentage': round(percentage_value, 3),
                    'min_max_normalised': min_max_normalised
                }

    return processed_dictionary


def process_alleles_and_groups(headers, rows, mhc_class):
    if mhc_class == 'class_i':
        fields = class_i_fields
    elif mhc_class == 'class_ii':
        fields = class_ii_fields
    else:
        fields = None

    if fields is None:
        return None

    allele_dict = {}
    allele_group_dict = {}

    super_populations = {}
    populations = {}

    sample_count = 0

    for row in rows:
        row_dict = dict(zip(headers, row))
        sample_id = row_dict[sample_id_field]
        super_population = row_dict[super_population_field]
        population = row_dict[population_field]
        
        update_simple_dictionary(super_populations, super_population)
        update_simple_dictionary(populations, population)

        sample_count += 1

        for field in row_dict:
            if field in fields:
                value = cleanup_allele_numbers(row_dict[field])
                if value is not None:
                    locus = field.split(' ')[0]
                    allele_group = f"{locus}*{value.split(':')[0]}"
                    allele_number = f"{locus}*{value}"

                    allele_group_dict = update_dictionary(allele_group_dict, allele_group, super_population, population, 'allele_group')
                    allele_dict = update_dictionary(allele_dict, allele_number, super_population, population, 'allele')


    processed_allele_groups = post_process_dictionary(allele_group_dict, populations, super_populations)
    processed_alleles = post_process_dictionary(allele_dict, populations, super_populations)
    
    return processed_alleles, processed_allele_groups, sample_count



                    
