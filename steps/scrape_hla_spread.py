import requests
from bs4 import BeautifulSoup

import os
import json

from common.helpers import deslugify_allele_group, slugify


header_mapping = {
    'PMID': 'pubmed_id', 
    'Allele': 'allele_group', 
    'Disease': 'disease', 
    'Population': 'population', 
    'Drug Names': 'drug_names', 
    'SNP': 'snp', 
    'Class': 'classification', 
    'Sentence': 'supporting_sentence'
}

def parse_html_data(html:str, allele_group:str):
    disease_associations = {}
    related_alleles = {}

    data = []
    soup = BeautifulSoup(html, 'html.parser')

    table = soup.find(id='hla_fullreport')

    if table:
        table_header = table.find('thead')
        header_row = table_header.find('tr')
        headers = [header.text for header in header_row.find_all('th')]
        
        mapped_headers = [header_mapping.get(header, header) for header in headers]

        table_body = table.find('tbody')
        rows = table_body.find_all('tr')
        print (len(rows))
        for row in rows:
            raw_row = dict(zip(mapped_headers, [cell.text.strip() for cell in row.find_all('td')]))

            if 'allele_group' in raw_row:
                if raw_row['allele_group'] == allele_group:
                    data.append(raw_row)
                else:
                    if '*' in raw_row['allele_group']:
                        locus = raw_row['allele_group'].split('*')[0]
                        if locus in ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G']:
                            if locus not in related_alleles:
                                related_alleles[locus] = {}
                            if raw_row['allele_group'] not in related_alleles[locus]:
                                related_alleles[locus][raw_row['allele_group']] = {'count': 0 }
                            related_alleles[locus][raw_row['allele_group']]['count'] += 1
    print (len(data))
    for row in data:
        disease = row['disease']
        if 'positive' in row['classification']:
            if disease not in ['NA', 'disease', 'diseases', 'disease_progression', 'inflammation', 'infection', 'strains', 'strain', 'death', 'deaths', 'wound', 'relapse', 'recurrence', 'rs1057141', 'rs492602', 'rs2287987', 'rs3823342', 'rs10484555', 'rs17179220', 'rs9260997']:
                disease_slug = slugify(disease)
                if disease_slug not in disease_associations:
                    disease_associations[disease_slug] = {'count': 0, 'rows': [], 'term': disease}
                disease_associations[disease_slug]['count'] += 1
                disease_associations[disease_slug]['rows'].append(row)

    return disease_associations, related_alleles


def fetch_data(allele_group_slug:str):
    url = f"https://hla-spread.igib.res.in/search/filter?limit=1000&search={allele_group_slug}"

    tmp_path = f"tmp/hla_spread/{allele_group_slug}.html"

    if not os.path.exists(tmp_path):
        response = requests.get(url)
        if response.status_code != 200:
            print(f"Failed to fetch {url}")
        else:
            os.makedirs(os.path.dirname(tmp_path), exist_ok=True)
            html = response.text
            with open(tmp_path, 'w') as f:
                f.write(html)
    else:
        with open(tmp_path, 'r') as f:
            html = f.read()

    return html


def get_disease_associations(allele_group_slug:str):
    html = fetch_data(allele_group_slug)

    allele_group = deslugify_allele_group(allele_group_slug)

    if html:
        disease_associations, related_alleles = parse_html_data(html, allele_group)

        # sort the related_alleles associations by count
        for locus in related_alleles:
            related_alleles[locus] = dict(sorted(related_alleles[locus].items(), key=lambda item: item[1]['count'], reverse=True))
        

        # sort the disease associations by count
        disease_associations = dict(sorted(disease_associations.items(), key=lambda item: item[1]['count'], reverse=True))

        for disease in disease_associations:
            print (f"{disease}: {disease_associations[disease]['count']}")

        with open(f"tmp/hla_spread/{allele_group_slug}.json", 'w') as f:
            f.write(json.dumps(disease_associations, indent=4))
    
    return disease_associations, related_alleles


def scrape_hla_spread(locus:str):
    print (f"Scraping HLA Spread for {locus}")

    locus_slug = locus.lower().replace('-', '_')


    input_filename = f"output/processed_data/allele_groups/{locus_slug}.json"

    with open(input_filename, 'r') as f:
        allele_groups = json.load(f)

    locus_associations = {}

    for allele_group_slug in allele_groups.keys():
        print (allele_group_slug)
        if allele_group_slug not in locus_associations:
            locus_associations[allele_group_slug] = {}
        
        disease_associations, related_alleles = get_disease_associations(allele_group_slug)

        locus_associations[allele_group_slug] = {'disease_associations': {key: value['count'] for key, value in disease_associations.items()}, 'total_count': len(disease_associations)}

    print (locus_associations)

    with open(f"output/processed_data/hla_spread/{locus_slug}_spread.json", 'w') as f:
        f.write(json.dumps(locus_associations, indent=4))

    return locus_associations


if __name__ == '__main__':
    hla_spread = {}
    for locus in ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G']:
        locus_slug = slugify(locus)
        hla_spread[locus_slug] = scrape_hla_spread(locus)
        
    with open("output/processed_data/hla_spread/hla_spread.json", 'w') as f:
        f.write(json.dumps(hla_spread, indent=4))
    


    
