
from common.pipeline import load_config

from fetch_ipd_data import fetch_ipd_data
from fetch_h2_data import fetch_h2_data
from parse_class_i_locus_data import construct_class_i_locus_allele_lists
from parse_class_i_bulk_data import construct_class_i_bulk_allele_lists
from create_folder_structure import create_folder_structure


from rich.console import Console
console = Console()

from rich import print


steps = [
    create_folder_structure,
    fetch_ipd_data,
    fetch_h2_data,
    construct_class_i_locus_allele_lists,
    construct_class_i_bulk_allele_lists,
    construct_class_i_locus_allele_lists
]

def step_number(step:str):
    return int(step) - 1

def run_step(*args, **kwargs):
    step = step_number(args[0])
    args = args[1:]
    log_data = steps[step](*args, **kwargs)
    print (log_data)
    print (args)
    print (kwargs)
    pass



action_log = {}

print ("")
console.rule(title="Running Allele Pipeline")
print ("")

console.print(f"There are {len(steps)} steps to this pipeline")
console.print("1. Creating the folder structure in the output directory")
console.print("2. Downloading latest versions of the IPD sequence sets (for human and most non-human species)")
console.print("3. Downloading latest versions of the mouse sequence sets")
console.print("4. Parsing the human Class I sequences from IPD")
console.print("5. Parsing the non-human Class I sequences from IPD")
console.print("6. Parsing the mouse Class I sequences")

print ("")
config = load_config()

console.rule(title="1. Creating the folder structure in the output directory")

create_folder_structure(config)

console.rule(title="2. Retrieving IPD sequence sets")

fetch_ipd_data(config)

console.rule(title="3. Retrieving mouse sequence sets")

fetch_h2_data(config)

hla_class_i = config['HLA_CLASS_I']

i = 1
for locus in hla_class_i:
    console.rule(title=f"4.{i}. Parsing IPD sequence set for HLA-{locus}")

    construct_class_i_locus_allele_lists(locus, 'hla', 'IPD_IMGT_HLA_PROT')

    i+=1

console.rule(title=f"5. Parsing IPD sequence set for non-human Class I")

construct_class_i_bulk_allele_lists('IPD-MHC')

h2_class_i = config['H2_CLASS_I']

j = 1
for locus in h2_class_i:
    console.rule(title=f"6.{j}. Parsing custom sequence set for H2-{locus}")

    construct_class_i_locus_allele_lists(locus, 'h2', 'H2_CLASS_I_PROT')

    j+=1

run_step("1",config,verbose=True)

