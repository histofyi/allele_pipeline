# Allele Pipeline

A pipeline to manipulate MHC Class I and II sequences from IPD resources

## Setup

The pipeline will need to have three directories created in the root of the folder:

- `tmp`
- `output`
- `logs`

These are not included in this repository as they will be filled with the outputs of steps. They are automatically created if they don't already exist in the first step of the pipeline.

The pipeline also needs some TOML formatted configuration files which are modular and assembled into a confg dictionary by the load_config() function in common/pipeline.py.

There are four files in total. 

### config.toml

Included in the repository. This file contains the filenames for the other configuration files.

### constants.toml

Included in the repository. This file contains constants used by the pipeline.

### paths.toml

Not included in the repository as it will depend on each user's directory structure.

Parameters for the following keys must be included.

- `WAREHOUSE_PATH`
- `PIPELINE_PATH`
- `TMP_PATH`
- `OUTPUT_PATH`
- `LOG_PATH`

### secrets.toml

Currently not yet used. Will be used for uploading logs and compiled files to AMAZON S3 for the histo.fyi implementation of this pipeline.