import os
import datetime

from typing import Dict
from rich import print, rule

import requests

from common.pipeline import load_config

def fetch_h2_data(config:Dict, force:bool = False):

    print("")
    for datasource in ['H2_CLASS_I_PROT']:
        
        filepath = f'{config["TMP_PATH"]}/{datasource.lower()}.fasta'

        if not os.path.exists(filepath) or force is True:

            fasta_url = config[datasource]

            print (f"Fetching {datasource} sequences")
            print ('')

            print (fasta_url)

            download_command = f'wget {fasta_url} -O {filepath}'
            os.system(download_command)

        else:

            print (f"{datasource} sequences already downloaded on {datetime.datetime.fromtimestamp(os.path.getmtime(filepath))}")
            print ('')


        
