import json
import os

import requests
from bs4 import BeautifulSoup
import country_converter as coco


def convert_slug_to_allele_frequencies(allele_slug):
    allele_components = allele_slug.split('_')
    return f"{allele_components[1]}*{allele_components[2]}:{allele_components[3]}".upper()


def fetch_allele_frequencies(allele_slug, num_rows=None):
    cached = False

    cache_location = f"tmp/allele_frequencies/{allele_slug}.json"

    if os.path.exists(cache_location):
        with open(cache_location, 'r') as f:
            rows = json.load(f)
            cached = True

    if not cached:
        
        url = f"https://www.allelefrequencies.net/hla6002a.asp?all_name={convert_slug_to_allele_frequencies(allele_slug)}"
        rows = []

        r = requests.get(url)

        if r.status_code == 200:
            html = r.content

            soup = BeautifulSoup(r.content, 'html.parser')

            allele_table_rows = soup.find('table', {'class': 'tblNormal'}).find_all('tr')

            for i, row in enumerate(allele_table_rows):
                if i != 0:
                    row_list = [td for td in (row.find_all('td'))]
                    if len(row_list[2]) > 0:
                        threeletter = row_list[0].find('img')['alt']
                        twoletter = coco.convert(names=threeletter, to='ISO2')
                        row_dict = {
                            'threeletter': threeletter,
                            'twoletter': twoletter,
                            'population': row_list[1].text,
                            'phenotype_frequency': float(row_list[2].text),
                            'allele_frequency': float(row_list[3].text),
                            'sample_size': int(row_list[5].text)
                        }
                        rows.append(row_dict)

            with open(cache_location, 'w') as f:
                json.dump(rows, f, indent=4)

    if num_rows:
        return rows[:num_rows]
    else:
        return rows


locus_slug = 'hla_g'

input_filename = f"output/processed_data/allele_groups/{locus_slug}.json"

allele_groups = json.load(open(input_filename, 'r'))

for allele_group in allele_groups:

    for allele_slug in allele_groups[allele_group]:
        print (allele_slug)

        rows = fetch_allele_frequencies(allele_slug, num_rows=20)


