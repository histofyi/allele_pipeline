from typing import Dict, List, Optional, Union

import toml
import json

import git
from importlib.metadata import version
from dparse import parse, filetypes
        
import datetime
import hashlib

from rich import print

def load_config(console, verbose:bool=False) -> Dict:
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
        console.print("Configuration")
        console.print (config)
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


def get_dependencies() -> Dict:
    with open('Pipfile','r') as file:
        pipfile = parse(file.read(), file_type=filetypes.pipfile)
    json_pipfile = json.loads(pipfile.json())
    dependencies = [dependency['name'] for dependency in json_pipfile['dependencies']]
    versions = [version(dependency) for dependency in dependencies]
    return {k:v for k,v in zip(dependencies,versions)}
        


class Pipeline():
    def __init__(self, steps:Dict, console, verbose:bool=False, mode:str='development', force:bool=False):

        self.logoutput = True
        self.verbose = verbose
        self.force = force
        self.console = console
        self.steps = steps
        self.mode = mode

        self.config = load_config(self.console, verbose=self.verbose )
        
        self.get_repository_info()
        self.initialise()

        
    def run_step(self, step_number, **kwargs):
        
        if 'substep' in kwargs:
            step_title_number = f"{step_number}.{kwargs['substep']}"
            substep = kwargs['substep']
        else:
            step_title_number = str(step_number)
            substep = None
        
        self.console.rule(title=f"{step_title_number}. {self.steps[step_number]['title_template'].format(**kwargs)}")
        
        started_at = get_current_time()

        kwargs['verbose'] = self.verbose
        kwargs['force'] = self.force
        kwargs['output_path'] = self.output_path
        kwargs['log_path'] = self.log_path

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
        

    def get_repository_info(self):
        repo = git.Repo(search_parent_directories=True)
        self.repository_name = repo.remotes.origin.url.split('.git')[0].split('/')[-1]
        self.pipeline_version = repo.head.object.hexsha
        self.pipeline_name = self.repository_name.replace('_',' ').capitalize()


    def initialise(self):
        if self.mode == 'release':
            # switch the output directory to the warehouse
            self.output_path = f"{self.config['WAREHOUSE_PATH']}/{self.config['PIPELINE_WAREHOUSE_FOLDER']}"
            self.log_path = f"{self.config['WAREHOUSE_PATH']}/logs/{self.config['PIPELINE_WAREHOUSE_FOLDER']}"
        else:
            self.output_path = self.config['OUTPUT_PATH']
            self.log_path = self.config['LOG_PATH']

        started_at = get_current_time()
        self.action_logs = {
            'started_at': started_at,
            'steps':{},
            'repository_name': self.repository_name,
            'pipeline_name': self.pipeline_name,
            'pipeline_version': self.pipeline_version
        }
        
        self.console.print ("")
        self.console.rule(title="Initialising...")
        self.console.print ("")
        self.console.print (f"{self.pipeline_name} (commit sha : {self.pipeline_version}) started at {started_at}")
        self.console.print ("")
        self.console.rule(title=f"Running {self.pipeline_name}")
        self.console.print ("")
        self.console.print(f"There are {len(self.steps)} steps to this pipeline")
        for step in self.steps:
            self.console.print(f"{step}. {self.steps[step]['list_item']}")
        self.console.print("")


    def finalise(self):
        self.action_logs['dependencies'] = get_dependencies()
        self.action_logs['completed_at'] = get_current_time()
        start = datetime.datetime.fromisoformat(self.action_logs['started_at'])
        end = datetime.datetime.fromisoformat(self.action_logs['completed_at'])
        delta = end - start
        datehash = hashlib.sha256(self.action_logs['completed_at'].encode('utf-8')).hexdigest()
        logfilename = f"{self.log_path}/{self.repository_name}-{datehash}.json"
        with open(logfilename, 'w') as logfile:
            logfile.write(json.dumps(self.action_logs, sort_keys=True, indent=4))
        self.console.print(f"Pipeline completed at {self.action_logs['completed_at']} : Execution time : {delta}") 
        return self.action_logs
