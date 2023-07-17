
from common.pipeline import load_config, Pipeline, run_step

from parse_class_i_locus_data import construct_class_i_locus_allele_lists
from parse_class_i_bulk_data import construct_class_i_bulk_allele_lists
from fetch_raw_data import fetch_raw_datasets
from create_folder_structure import create_folder_structure


from rich.console import Console
console = Console()

from rich import print

import git

steps = {
    '1':{
        'function':create_folder_structure,
        'title_template':'Creating the folder structure in the output directory',
        'list_item':'Creating the folder structure in the output directory'
    },
    '2':{
        'function':fetch_raw_datasets,
        'title_template':'Downloading latest versions of the IPD and H2 sequence datasets',
        'list_item':'Downloading latest versions of the IPD and H2 sequence datasets'
    },
    '3':{
        'function':construct_class_i_locus_allele_lists,
        'title_template':'{substep}. Parsing IPD sequence set for HLA-{locus}',
        'list_item':'Parsing the human Class I sequences from IPD'
    },
    '4':{
        'function':construct_class_i_bulk_allele_lists,
        'title_template':'Parsing IPD sequence set for non-human Class I',
        'list_item':'Parsing the non-human Class I sequences from IPD'
    },
    '5':{
        'function':construct_class_i_locus_allele_lists,
        'title_template':'',
        'list_item':'{substep}. Parsing H2 sequence set for H2-{locus}'
    }
}





pipeline = Pipeline(steps, console, logoutput=True)

config = load_config()

pipeline.run_step('1')

pipeline.run_step('2')

hla_class_i = config['HLA_CLASS_I']

i = 1
for locus in hla_class_i:
    pipeline.run_step('3', substep=i, locus=locus, species_slug='hla', sequence_set='IPD_IMGT_HLA_PROT')
    i+=1

console.rule(title=f"4. ")

construct_class_i_bulk_allele_lists('IPD-MHC')

h2_class_i = config['H2_CLASS_I']

j = 1
for locus in h2_class_i:
    pipeline.run_step('5', substep=j, locus=locus, species_slug='h2', sequence_set='H2_CLASS_I_PROT')
    j+=1



