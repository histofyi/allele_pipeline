from typing import Dict

import os
from pathlib import Path

from rich.console import Console
console = Console()


def create_folder(folder_path:str, verbose:bool) -> str:
    """
    This function creates the folder path if it does not already exist 
    
    Args:
        folder_path (str): the full path to the folder
        verbose (bool): whether the function should echo to the terminal if this argument is set to True
    
    Returns:
        str: the status of the folder e.g. created, already in existence
    """

    # check if the folder exists
    if not os.path.exists(folder_path):
        # if it doesn't exist, set folder_status to `folders_created`
        folder_status = 'folders_created'
        # create the folder and any parent folders needed
        Path(folder_path).mkdir(parents=True, exist_ok=True)
        # if verbose is set to True, send a message to the terminal
        if verbose:  
            console.print (f"{folder_path} created")  
    else:
        # if it does exist, set folder_status to `folders_in_existence`
        folder_status = 'folders_in_existence'

        # if verbose is set to True, send a message to the terminal
        if verbose:
            console.print (f"{folder_path} already exists")  
    return folder_status


def create_folder_structure(config, **kwargs) -> Dict:
    """
    This function creates the folder structure for the outputs of the pipeline

    Args:
        config (Dict): the configuration from the config.toml file
        verbose (bool): whether this step should output to the terminal or not (in kwargs)

    Returns:
        Dict : a dictionary of actions performed 
    """

    # extract appropriate variables from kwargs
    verbose = kwargs['verbose']
    output_path = kwargs['output_path']
    log_path = kwargs['log_path']

    # default action log file for this step
    action_log = {
        'folders_created':[],
        'folders_in_existence':[], 
        'completed_at': None
    }

    # create a list of folders from the config
    folders = [folder for folder in config['SEQUENCE_TYPES']]
    # append the `protein_alleles` folder
    folders.append('protein_alleles')

    # iterate through the folders in the list
    for folder in folders:
        folder_path = f"{output_path}/processed_data/{folder}"
        folder_status = create_folder(folder_path, verbose)
        action_log[folder_status].append(folder)

    # create the local logfile folder
    folder_status = create_folder(log_path, verbose)
    action_log[folder_status].append('log')

    # create the local tmp folder
    folder_status = create_folder(config['TMP_PATH'], verbose)
    action_log[folder_status].append('tmp')

    return action_log



