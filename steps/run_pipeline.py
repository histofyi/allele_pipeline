from typing import Dict

from common.pipeline import Pipeline

from parse_class_i_locus_data import construct_class_i_locus_allele_lists
from parse_class_i_bulk_data import construct_class_i_bulk_allele_lists
from fetch_raw_data import fetch_raw_datasets
from create_folder_structure import create_folder_structure

from rich.console import Console

def run_pipeline(verbose:bool=False, force:bool=False, logoutput:bool=True) -> Dict:
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
            'title_template':'Parsing IPD sequence set for HLA-{locus}',
            'list_item':'Parsing the human Class I sequences from IPD'
        },
        '4':{
            'function':construct_class_i_bulk_allele_lists,
            'title_template':'Parsing IPD sequence set for non-human Class I',
            'list_item':'Parsing the non-human Class I sequences from IPD'
        },
        '5':{
            'function':construct_class_i_locus_allele_lists,
            'title_template':'Parsing H2 sequence set for H2-{locus}',
            'list_item':'Parsing the mouse Class I sequences from a custom dataset'
        }
    }

    pipeline = Pipeline(steps, Console(), logoutput)

    pipeline.run_step('1')

    pipeline.run_step('2')

    hla_class_i = pipeline.get_config_item('HLA_CLASS_I')

    i = 1
    for locus in hla_class_i:
        pipeline.run_step('3', substep=i, locus=locus, species_slug='hla', sequence_set='IPD_IMGT_HLA_PROT')
        i+=1

    construct_class_i_bulk_allele_lists('IPD-MHC')

    h2_class_i = pipeline.get_config_item('H2_CLASS_I')

    j = 1
    for locus in h2_class_i:
        pipeline.run_step('5', substep=j, locus=locus, species_slug='h2', sequence_set='H2_CLASS_I_PROT')
        j+=1
    
    return {}


def main():
    output = run_pipeline()

main()



