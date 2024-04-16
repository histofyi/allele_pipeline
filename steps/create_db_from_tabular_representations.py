from typing import Dict, List, Tuple

import os
import json
import csv

from common.helpers import write_csv_file

def create_db_from_csvs(input_filenames:List[str], output_filename:str) -> None:
    if os.path.exists(output_filename):
        os.remove(output_filename)

    if len(input_filenames) > 1:
        input_filename_str = ' '.join(input_filenames)
    else:
        input_filename_str = input_filenames[0]

    csvs_to_sqllite_command = f"csvs-to-sqlite {input_filename_str} {output_filename}"
    os.system(csvs_to_sqllite_command)

def combine_csv_files(filenames:List) -> List:
    i = 0
    table = []
    for filename in filenames:
        print (f"Appending {filename} to table")
        with open(filename, 'r') as filehandle:
            reader = csv.reader(filehandle)
            if i == 0:
                table = [row for row in reader]
                headers = table[0]
                data_row_length = len(table) - 1
                print (f"Appending headers {headers} and {data_row_length} data rows")
            data_rows = [row for row in reader][1:]
            print (f"Appending {len(data_rows)} data rows...")
            table = table + data_rows
            i += 1
        print (f"New table length {len(table)}")
    print (f"Completed combining files, table length {len(table)}")    
    return table






def create_db_from_tabular_representations(config:Dict, **kwargs) -> None:
    """
    This function takes a set of loci and builds a database table from them.

    Args:
        loci (List): the list of loci to be processed
        verbose (bool): whether specific information is output to the terminal, for large sequence sets this can be overwhelming and significantly slow down the function
    Returns:

    """
    loci = kwargs['loci']
    

    # combine allele sequence information

    allele_sequence_filenames = []

    for locus in loci:
        allele_sequence_filenames.append(f"output/tabular_data/alleles/{locus.lower()}.csv")
    table = combine_csv_files(allele_sequence_filenames)

    alleles_output_filename = "output/tabular_data/alleles.csv"

    write_csv_file(alleles_output_filename , table)

    create_db_from_csvs([alleles_output_filename], "output/tabular_data/alleles.db")


    allele_relationship_filenames = []

    for mode in ['motif','structure']:
        for locus in loci:
            if locus not in ['hla_e', 'hla_f', 'hla_g']:
                allele_relationship_filenames.append(f"output/tabular_data/relationships/{locus.lower()}_{mode}.csv")
    
    table = combine_csv_files(allele_relationship_filenames)

    relationships_output_filename = "output/tabular_data/relationships.csv"

    write_csv_file(relationships_output_filename , table)

    create_db_from_csvs([relationships_output_filename], "output/tabular_data/relationships.db")


    create_db_from_csvs([alleles_output_filename, relationships_output_filename], "output/tabular_data/combined.db")

    pass
