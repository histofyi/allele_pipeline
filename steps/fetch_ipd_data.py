import os
from rich import print, rule

from common.pipeline import load_config


config = load_config()

print ()

for datasource in ['IPD_IMGT_HLA_PROT','IPD_MHC_PROT']:

    
    filepath = f'{config["TMP_PATH"]}/{datasource.lower()}.fasta'
    fasta_url = config[datasource]

    print (datasource)
    print (filepath)
    print (fasta_url)
    print ('')
    rule

    download_command = f'wget {fasta_url} -O {filepath}'

    os.system(download_command)

    #os.system()