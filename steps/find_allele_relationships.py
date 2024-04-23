from typing import Dict, List, Tuple

from Levenshtein import hamming

import json
import csv

def locate_polymorphisms(allele_pseudosequence:str, match_pseudosequence:str, pocket_positions:List) -> Dict:
    polymorphisms = {}
    for i, position in enumerate(pocket_positions):
        if allele_pseudosequence[i] != match_pseudosequence[i]:
            polymorphisms[position] = {'from':match_pseudosequence[i], 'to':allele_pseudosequence[i]}
    return polymorphisms

def find_closest_alleles_function(allele_slug:str, known_alleles:List[str], known_pseudosequences:List, allele_pseudosequence:str, pocket_positions:List, mode:str='motif') -> Dict:
    closest_alleles = None
    
    # the experimental label is set to true if the allele slug matches that of a known motif or structure. 
    # the nearest label is set to the allele_slug of the nearest allele which has a known motif or structure
    nearest_label = f"nearest_known_allele"

    # first we check if the allele_slug matches one of those in the known alleles
    if allele_slug in known_alleles:
        closest_alleles = [{
            'nearest_known_allele':allele_slug,
            'distance':0,
            'relationship_label':f"experimentally_determined_{mode}",
            'polymorphisms': None
        }]
    # if it doesn't, we'll first check if the allele pseudosequence 
    else:
        
        # we set the initial minimum score to the length of the pseudosequence, we're using the hamming distance so the initial value would be a sequence with no matches at all
        min_score = len(allele_pseudosequence)
        # we create a lowest scores dictionary, we want to collect matching alleles with the lowest scores, but we won't know this as we iterate through the known alleles
        lowest_scores = {}

        # we iterate through the known alleles
        for known_allele in known_pseudosequences:
            # if the pseudosequence matches the known pseudosequence exactly, we can stop the search and return the known allele
            if allele_pseudosequence == known_pseudosequences[known_allele]:
                closest_alleles = [{
                    'nearest_known_allele':known_allele,
                    'distance':0,
                    'relationship_label':'exact_pseudosequence_match',
                    'polymorphisms': None                
                }]
                break
            # if the pseudosequence doesn't match exactly, we calculate the hamming distance between the pseudosequence for the allele and for the known test allele
            else:
                score = hamming(allele_pseudosequence, known_pseudosequences[known_allele])
                # if the score is less than or equal to the minimum score, we update the minimum score and add the allele to the lowest scores dictionary
                if score <= min_score or score in lowest_scores:
                    polymorphisms = locate_polymorphisms(allele_pseudosequence, known_pseudosequences[known_allele], pocket_positions)
                    if score < min_score:
                        min_score = score
                    # if the score isn't already in the lowest scores dictionary, we add it
                    if score not in lowest_scores:
                        lowest_scores[score] = []
                    # we append the allele to the appropriate portion of the lowest scores dictionary
                    lowest_scores[score].append({
                        'nearest_known_allele':known_allele,
                        'distance':score,
                        'relationship_label':"nearest_pseudosequence_match",
                        'polymorphisms': polymorphisms
                    })
        if not closest_alleles:
            closest_alleles = lowest_scores[min_score]
        
    return closest_alleles


def tabulate_relationships(relationships:Dict, relationship_type:str, pocket_positions:List, distance_frequency_cutoff:int=10) -> List:
    rows = []
    outliers = []

    

    labels = ['allele_slug', 'known_allele_slug', 'distance', 'relationship_label', 'relationship_type']

    for position in pocket_positions:
        labels.append(f"alpha_{position}")


    rows.append(labels)
    
    for allele in relationships:
        if relationships[allele][0]['distance'] > distance_frequency_cutoff:
            outliers.append(allele)
        else:
            for relationship in relationships[allele]:
                row = [allele, relationship["nearest_known_allele"], relationship['distance'], relationship['relationship_label'], relationship_type]
                for position in pocket_positions:
                    if relationship['polymorphisms'] and position in relationship['polymorphisms']:
                        row.append(f"{relationship['polymorphisms'][position]['from']}{position}{relationship['polymorphisms'][position]['to']}")
                    else:
                        row.append('')
                rows.append(row)
    return rows, outliers


def find_allele_relationships(config:Dict, **kwargs) -> None:
    """
    This function takes the NetMHCPan sequence of each allele and checks for the Hamming distance to the closest allele

    Args:
        locus (str): the locus to be parsed
        verbose (bool): whether specific information is output to the terminal, for large sequence sets this can be overwhelming and significantly slow down the function
    Returns:

    """
    test_locus = kwargs['locus']
    loci = kwargs['loci']

    relationship_types = config['CONSTANTS']['RELATIONSHIP_TYPES']

    known_alleles = {}

    for relationship_type in relationship_types:
        config_key = f"{relationship_type.upper()}_ALLELES"
        known_alleles[relationship_type] = config['CONSTANTS'][config_key]

    pocket_positions = config['CONSTANTS']['IMGT_POCKET_RESIDUES']

    # we'll initialise some datastructures that we'll use to store the pseudosequences for the known motifs and structures and the alleles we want to test

    alleles_to_test = None

    pseudosequences = {}

    for relationship_type in relationship_types:
        pseudosequences[relationship_type] = {}

    # we'll iterate through the loci 
    for locus in loci:
        # we'll load the pseudosequences for the known motifs and structures
        input_filename = f"output/processed_data/protein_alleles/{locus.lower()}.json"
        with open(input_filename, 'r') as protein_alleles_filehandle:
            raw_alleles = json.load(protein_alleles_filehandle)

        # we'll iterate through the alleles in the raw alleles
        for allele in raw_alleles:
            
            for relationship_type in relationship_types:
                # we'll check if the allele is in the known motifs
                if allele in known_alleles[relationship_type]:
                    pseudosequences[relationship_type][allele] = raw_alleles[allele]['pocket_pseudosequence']
                    
       # if the locus is the locus we're testing, we'll add the alleles to the alleles to test dictionary
        if locus == test_locus:
            alleles_to_test = {}

            # we'll iterate through the alleles in the raw alleles
            for allele in raw_alleles:
                canonical_protein_allele_name = raw_alleles[allele]['canonical_allele']['protein_allele_name']
                # we'll check if the canonical protein allele name doesn't end in N or Q (these have differential or no expression and often contain deletions)
                if not canonical_protein_allele_name[-1] in ['N','Q']:
                    alleles_to_test[allele] = raw_alleles[allele]['pocket_pseudosequence']

    # we'll initialise some datastructures to store the related alleles



    print (f"Testing alleles for {test_locus}")

    related_alleles = {}
    outlier_alleles = {}

    for relationship_type in relationship_types:
        related_alleles[relationship_type] = {}
        outlier_alleles[relationship_type] = {}

    for allele_slug in alleles_to_test:
        for relationship_type in relationship_types:            
            related_alleles[relationship_type][allele_slug] = find_closest_alleles_function(allele_slug, known_alleles[relationship_type], pseudosequences[relationship_type], alleles_to_test[allele_slug],pocket_positions,  mode=relationship_type)

    
    # metrics needed

    # which alleles have the highest distance from the nearest known allele both structurally and from a motif standpoint
    # what is the distribution of distances
    # which known alleles are the most common nearest alleles
    # which alleles are related to the known alleles e.g. look up which alleles are identical to HLA-A*02:01 in NetMHCPan pseudosequence
    # which known alleles are related to the alleles e.g look up which known alleles are related to HLA-A*02:97


    # allele_slug, known_allele_slug, distance, relationship_type, mode


    # we'll set the distance frequency cutoff to 10 to check for outliers
    distance_frequency_cutoff = 10
    outlier_alleles = {}

    for mode in relationship_types:
        output_filename = f"output/tabular_data/relationships/{test_locus.lower()}_{mode}.csv"
        print (output_filename)
        table, outlier_alleles['motif'] = tabulate_relationships(related_alleles[mode], mode, pocket_positions, distance_frequency_cutoff)
        print (len(table))

        with open(output_filename, 'w', newline='\n') as f:
            writer = csv.writer(f)
            writer.writerows(table)

    print (f"Alleles with potential issues to explore: {outlier_alleles}")

    pass
