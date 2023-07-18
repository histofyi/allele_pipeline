from typing import Dict

import os
import datetime


from rich.console import Console
console = Console()

from common.helpers import slugify
from common.pipeline import get_current_time

def fetch_raw_datasets(config:Dict, force:bool=False, verbose:bool=False):

    action_log = {}
    for datasource in ['IPD_IMGT_HLA_PROT','IPD_MHC_PROT', 'H2_CLASS_I_PROT']:
        
        filepath = f'{config["TMP_PATH"]}/{datasource.lower()}.fasta'

        downloaded_at = None
        
        # TODO force a reload if the downloaded at date is longer ago than x days

        # if the file doesn't exist, or if we're forcing reload of sequences
        if not os.path.exists(filepath) or force is True:

            # pull the url from the config
            fasta_url = config[datasource]

            # if verbose is true, output to the terminal
            if verbose:
                print("")
                print (f"Fetching {datasource} sequences")
                print("")
            
            # downloard the file
            download_command = f'wget {fasta_url} -O {filepath}'
            os.system(download_command)

            action_log[slugify(datasource)] = {'downloaded':get_current_time()}

        else:
            
            downloaded_at = datetime.datetime.fromtimestamp(os.path.getmtime(filepath)).isoformat()

            if verbose:
                print("")
                print (f"{datasource} sequences already downloaded on {downloaded_at}")
                print("")
            
            action_log[slugify(datasource)] = {'cached':downloaded_at}
    
    return action_log



        
