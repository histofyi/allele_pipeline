# Allele Pipeline

A pipeline to manipulate MHC Class I and II sequences from IPD resources

## Setup

The pipeline will need to have three directories created in the root of the folder:

- `tmp`
- `output`
- `logs`

These are not included in this repository as they will be filled with the outputs of steps. 

The pipeline also needs a configuration file created `config.toml`

This must include the following parameters:

IPD_IMGT_HLA_PROT = "https://github.com/ANHIG/IMGTHLA/raw/Latest/hla_prot.fasta"
IPD_MHC_PROT = "https://github.com/ANHIG/IPDMHC/raw/Latest/MHC_prot.fasta"
H2_CLASS_I_PROT = "https://raw.githubusercontent.com/histofyi/datasets/main/h2_prot.fasta"

HLA_CLASS_I = ["A","B","C","E","F","G"]
H2_CLASS_I = ["K","D","L"]

WAREHOUSE_PATH = "/Users/username/code/warehouse"
PIPELINE_PATH = "/Users/username/code/allele_pipeline"
TMP_PATH = "/Users/username/code/allele_pipeline/tmp"

OUTPUT_PATH = "/Users/username/code/allele_pipeline/output"
LOG_PATH = "/Users/username/code/allele_pipeline/logs"

