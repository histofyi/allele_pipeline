from typing import Dict, Optional
import toml

import datetime
from rich import print

def load_config() -> Dict:
    """ 
    Loads the configuration file for the pipline and returns a dictionary of values

    Returns:
        dict: a dictionary of configuration variables 

    """
    config = toml.load('config.toml')
    return config


def read_lockfile(filename:str) -> Optional[str]:
    """
    Reads the content of a lockfile
    
    Args:
        filename (str) - the filename of the lockfile

    Returns:
        str: the contents of the 
    """
    try:
        with open(filename) as lockfile:
            test = lockfile.read()
    except:
        write_lockfile(filename, 'intialised')
        test = 'initialised'
    return test


def write_lockfile(filename:str, text:str):
    """
    Writes the content of a lockfile

    Args:
        filename (str) - the filename of the lockfile
        text (str) - the contents for the lockfile
    """
    with open(filename, 'w') as lockfile:
        lockfile.write(text)
    pass


def generate_completed_at_entry(action_log:Dict) -> Dict:
    # generate the completed at variable
    completed_at = datetime.datetime.now().isoformat()

    action_log['completed_at'] = completed_at

    # output completed at message to terminal
    print (f"Completed at {completed_at}")

    return action_log