from typing import List

import csv

def slugify(string:str) -> str:
    slug_char = '_'
    to_replace = [' ','-','.',',','[',']','{','}','(',')','/','\\','*',':']
    for replace_me in to_replace:
        if replace_me in string:
            string = string.replace(replace_me, slug_char)
    return string.lower()


def write_csv_file(filename:str, table:List):
    with open(filename, 'w', newline='\n') as filehandle:
        writer = csv.writer(filehandle)
        writer.writerows(table)  