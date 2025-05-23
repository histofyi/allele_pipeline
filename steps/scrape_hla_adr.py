import requests
from bs4 import BeautifulSoup

from common.helpers import slugify

import json
import os

tmp_dir = 'tmp/hla_adr'
filename = f"{tmp_dir}/results.html"
url = "https://www.allelefrequencies.net/hla-adr/adr_query.asp?dis_gene=&dis_allele=&dis_non_standard=&dis_drug=&dis_ethnicity=&dis_pvalue=&dis_disease=&dis_adr=&dis_country=&dis_geog_region=&dis_sort=&dummy=dummy"

adr_url_stem = 'https://www.allelefrequencies.net/hla-adr/'

if not os.path.exists(tmp_dir):
    os.makedirs(tmp_dir)


if not os.path.exists(filename):

    response = requests.get(url)

    if response.status_code == 200:
        html = response.text
        with open(filename, 'w') as f:
            f.write(html)

else:
    with open(filename, 'r') as f:
        html = f.read()

soup = BeautifulSoup(html, 'html.parser')

table = soup.find('table', {'class': 'hla_adr.tblNormal'})

rows = table.find_all('tr')

adverse_drug_reactions = {}

for row in rows:
    cells = [cell.text for cell in row.find_all('td')]
    links = row.find_all('a')
    if len(cells) > 0:
        pubmed_id = cells[1]
        drug = cells[2]
        allele = cells[3]
        ancestry = cells[5]
        if '*' in allele:
            locus_letter = allele.split('*')[0]
            if locus_letter in ['A', 'B', 'C', 'E', 'F', 'G']:
                locus = f"HLA-{locus_letter}"
                locus_slug = slugify(locus)

                allele_group = f"HLA-{allele.split(':')[0]}"
                allele_group_slug = slugify(allele_group)
                allele_number = f"HLA-{allele}"
                allele_slug = slugify(allele_number)
                print(pubmed_id, drug, locus, allele_number, allele_group, ancestry)
                for link in links:
                    if 'adr_report' in link['href']:
                        adr_url = adr_url_stem + link['href']


                if locus_slug not in adverse_drug_reactions:
                    adverse_drug_reactions[locus_slug] = {'locus': locus, 'allele_groups': {}, 'reaction_count': 0}
                adverse_drug_reactions[locus_slug]['reaction_count'] += 1
                
                if allele_group_slug not in adverse_drug_reactions[locus_slug]['allele_groups']:
                    adverse_drug_reactions[locus_slug]['allele_groups'][allele_group_slug] = {'allele_group': allele_group, 'alleles': {}, 'reaction_count': 0}
                adverse_drug_reactions[locus_slug]['allele_groups'][allele_group_slug]['reaction_count'] += 1

                if allele_slug not in adverse_drug_reactions[locus_slug]['allele_groups'][allele_group_slug]['alleles']:
                    adverse_drug_reactions[locus_slug]['allele_groups'][allele_group_slug]['alleles'][allele_slug] = {'allele_number': allele_number, 'reactions': [], 'reaction_count': 0}
                adverse_drug_reactions[locus_slug]['allele_groups'][allele_group_slug]['alleles'][allele_slug]['reactions'].append({
                        'pubmed_id': pubmed_id,
                        'drug': drug,
                        'ancestry': ancestry,
                        'adr_url': adr_url
                    })
                adverse_drug_reactions[locus_slug]['allele_groups'][allele_group_slug]['alleles'][allele_slug]['reaction_count'] += 1
            else:
                locus = None
            

print (adverse_drug_reactions)

with open('output/processed_data/hla_adr/hla_adr.json', 'w') as f:
    f.write(json.dumps(adverse_drug_reactions, indent=4))