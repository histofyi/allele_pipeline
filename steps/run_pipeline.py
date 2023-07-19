from typing import Dict

from common.pipeline import Pipeline

from parse_class_i_locus_data import construct_class_i_locus_allele_lists
from parse_class_i_bulk_data import construct_class_i_bulk_allele_lists
from fetch_raw_data import fetch_raw_datasets
from create_folder_structure import create_folder_structure

from rich.console import Console
import argparse

def run_pipeline(verbose:bool=False, force:bool=False, mode:str='development') -> Dict:
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

    pipeline = Pipeline(steps, Console(), force=force, verbose=verbose, mode=mode)

    pipeline.run_step('1')

    pipeline.run_step('2')

    hla_class_i = pipeline.get_config_item('HLA_CLASS_I')

    i = 1
    for locus in hla_class_i:
        pipeline.run_step('3', substep=i, locus=locus, species_slug='hla', sequence_set='IPD_IMGT_HLA_PROT')
        i+=1

    pipeline.run_step('4', sequence_set='IPD_MHC_PROT')
        
    h2_class_i = pipeline.get_config_item('H2_CLASS_I')

    j = 1
    for locus in h2_class_i:
        pipeline.run_step('5', substep=j, locus=locus, species_slug='h2', sequence_set='H2_CLASS_I_PROT')
        j+=1

    action_logs = pipeline.finalise()

    print (mode)
    
    return action_logs


def main():
    parser = argparse.ArgumentParser(prog='Allele Pipeline',
                    description='This pipeline builds a dataset of sequences, pseudosequences and canonical alleles.',
                    epilog='For more information see - https://github.com/histofyi/allele_pipeline')
    parser.add_argument('-v','--verbose', help='increases output verbosity (non-verbosity is the default)', action='store_true')
    parser.add_argument('-f', '--force', help='forces reloading of underlying datasets (not forcing reload is the default)', action='store_true')
    parser.add_argument('-r', '--release', help='switch between development and release modes (development mode is the default)', action='store_true')
    args = parser.parse_args() 

    if args.verbose:
        verbose = True
    else:
        verbose = False
    if args.force:
        force = True
    else:
        force = False
    if args.release:
        mode = 'release' 
    else:
        mode = 'development'    

    print (args)

    print (verbose)
    print (force)
    print (mode)
    output = run_pipeline(verbose=verbose, force=force, mode=mode)


if __name__ == '__main__':
    main()



