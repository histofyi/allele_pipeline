from typing import Dict, List, Optional, Union
import toml
import git

import datetime
from rich import print

def load_config(verbose:str=False) -> Dict:
    """ 
    Loads the configuration file for the pipline and returns a dictionary of values

    Returns:
        dict: a dictionary of configuration variables 

    """
    config = {}
    files = toml.load('config.toml')
    for file in files:
        this_config = toml.load(f"{files[file]}")
        for k,v in this_config.items():
            config[k] = v
    if verbose:
        print (config)
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


def get_current_time() -> str:
    return datetime.datetime.now().isoformat()



class Pipeline():
    def __init__(self, steps:Dict, console, verbose:bool=False, logoutput:bool=False, devmode:bool=True):
        repo = git.Repo(search_parent_directories=True)
        started_at = get_current_time()
        repository_name = repo.remotes.origin.url.split('.git')[0].split('/')[-1]
        pipeline_version = repo.head.object.hexsha
        pipeline_name = repository_name.replace('_',' ').capitalize()
        self.steps = steps
        self.config = load_config(verbose)
        self.action_logs = {
            'started_at': started_at,
            'steps':{},
            'repository_name': repository_name,
            'pipeline_name': pipeline_name,
            'pipeline_version': pipeline_version
        }
        self.console = console
        self.console.print ("")
        self.console.rule(title="Initialising...")
        self.console.print ("")
        self.console.print (f"{pipeline_name} (commit sha : {pipeline_version}) started at {started_at}")
        self.console.print ("")
        self.console.rule(title=f"Running {pipeline_name}")
        self.console.print ("")
        self.console.print(f"There are {len(self.steps)} steps to this pipeline")
        for step in steps:
            self.console.print(f"{step}. {self.steps[step]['list_item']}")
        self.console.print("")
        self.logoutput = logoutput


    def run_step(self, step_number, **kwargs):
        
        if 'substep' in kwargs:
            step_title_number = f"{step_number}.{kwargs['substep']}"
            substep = kwargs['substep']
        else:
            step_title_number = str(step_number)
            substep = None
        
        self.console.rule(title=f"{step_title_number}. {self.steps[step_number]['title_template'].format(**kwargs)}")
        
        started_at = get_current_time()

        _action_log = self.steps[step_number]['function'](self.config, **kwargs)

        if 'verbose' in kwargs:
            del kwargs['verbose']
        
        if 'substep' in kwargs:
            del kwargs['substep']

        completed_at = get_current_time()
        print (f"Completed at {completed_at}")

        self.action_logs['steps'][step_title_number] = {
            'step': step_number,
            'substep': substep,
            'step_action': self.steps[step_number]['function'].__name__,
            'started_at':started_at,
            'completed_at': completed_at,
            'arguments': kwargs,
            'action_log': _action_log
        }

        if self.logoutput:
            self.console.print(self.action_logs['steps'][step_title_number])
        pass

    def get_config_item(self, key:str) -> Union[None, str, int, List]:
        if key in self.config:
            return self.config[key]
        else: 
            return None
