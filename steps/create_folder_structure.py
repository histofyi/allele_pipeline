from typing import Dict

import os
from pathlib import Path

from rich.console import Console
console = Console()

from rich import print


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
            print (f"{folder_path} created")  
    else:
        # if it does exist, set folder_status to `folders_in_existence`
        folder_status = 'folders_in_existence'
        # if verbose is set to True, send a message to the terminal
        if verbose:
            print (f"{folder_path} already exists")  
    return folder_status


def create_folder_structure(config:Dict, verbose:bool=False) -> Dict:
    """
    This function creates the folder structure for the outputs of the pipeline

    Args:
        config (Dict): the configuration from the config.toml file
        verbose (bool): whether this step should output to the terminal or not

    Returns:
        Dict : a dictionary of actions performed 
    """

    # create a list of folders from the config
    folders = [folder for folder in config['SEQUENCE_TYPES']]
    
    # append the `protein_alleles` folder
    folders.append('protein_alleles')

    # default log file
    action_log = {
        'folders_created':[],
        'folders_in_existence':[]
    }

    # iterate through the folders
    for folder in folders:
        folder_path = f"{config['OUTPUT_PATH']}/{folder}"
        folder_status = create_folder(folder_path, verbose)
        action_log[folder_status].append(folder)


    folder_status = create_folder(config['LOG_PATH'], verbose)
    action_log[folder_status].append('log')

    folder_status = create_folder(config['TMP_PATH'], verbose)
    action_log[folder_status].append('tmp')

    return action_log



