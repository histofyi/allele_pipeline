
from common.pipeline import load_config

from fetch_ipd_data import fetch_ipd_data
from parse_class_i_imgt_data import construct_class_i_hla_allele_lists

from rich.console import Console
console = Console()

from rich import print

print ("")
console.rule(title="Running Allele Pipeline")
print ("")

console.print("There are 5 steps to this pipeline")
console.print("1. Downloading latest versions of the IPD sequence sets (for human and most non-human species)")
console.print("2. Downloading latest versions of the mouse sequence sets")

console.print("3. Parsing the human Class I sequences from IPD")
console.print("4. Parsing the non-human Class I sequences from IPD")
console.print("5. Parsing the mouse Class I sequences")

print ("")
config = load_config()

console.rule(title="1. Retrieving IPD sequence sets")

fetch_ipd_data(config)

console.rule(title="2. Retrieving mouse sequence sets - coming soon")

print ("")

hla_class_i = config['HLA_CLASS_I']

i = 1
for locus in hla_class_i:
    console.rule(title=f"3.{i}. Parsing IPD sequence set for HLA-{locus}")

    construct_class_i_hla_allele_lists(locus)

    i+=1



