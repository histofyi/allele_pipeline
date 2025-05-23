from typing import Dict, List, Tuple

import base64
from io import BytesIO

from common.helpers import slugify

from matplotlib.figure import Figure
import numpy as np

import json




def generate_allele_groups(pseudosequences:Dict) -> Dict:
    """
    This function takes a dictionary of pseudosequences and returns a dictionary of allele groups containing those pseudosequences.

    Args:
        pseudosequences (Dict): A dictionary of pseudosequences.
    
    Returns:
        Dict: A dictionary of allele groups containing pseudosequences.
    """
    #NOTE - we can't assume that the pseudosequences map only to a single allele group as some of them clearly do not

    # initialize the allele groups dictionary
    allele_groups = {}

    # iterate over the pseudosequences
    for pseudosequence in pseudosequences:

        # for each pseudosequence, we're going to generate a list of protein alleles associated with it
        # we don't care about the gene alleles, just the protein alleles, so we'll filter out the gene alleles
        protein_alleles = list(set([allele['protein_allele_name'] for allele in pseudosequences[pseudosequence]['alleles']]))

        # iterate over the unique protein alleles
        for allele in protein_alleles:
            # we're going to split the allele name on the colon to get the allele group
            allele_group = allele.split(':')[0]
            allele_group_slug = slugify(allele_group)

            # if the allele group isn't in the allele groups dictionary, we'll add it
            if allele_group_slug not in allele_groups:
                allele_groups[allele_group_slug] = {
                    'allele_group':allele_group,
                    'pseudosequences':{}
                }

            # if the pseudosequence isn't in the allele group, we'll add it    
            if pseudosequence not in allele_groups[allele_group_slug]['pseudosequences']:
                allele_groups[allele_group_slug]['pseudosequences'][pseudosequence] = {'count':0, 'alleles':[]}
            # and we'll add the allele to the pseudosequence and increment the count
            allele_groups[allele_group_slug]['pseudosequences'][pseudosequence]['alleles'].append(allele)
            allele_groups[allele_group_slug]['pseudosequences'][pseudosequence]['count'] += 1

    return allele_groups


def find_lowest_allele_number(alleles:List) -> str:
    """
    This function takes a list of alleles and returns the allele with the lowest allele number

    Args:
        alleles (List): A list of alleles
    
    Returns:
        str: The allele number of the allele with the lowest allele number in the list
    """
    # we're going to split the allele on the colon and get the allele number
    allele_number_numerics = [int(allele.split(':')[1]) for allele in alleles]
    # we're going to split the allele on the colon and get the allele group
    allele_group = alleles[0].split(':')[0]
    # we're going to find the lowest allele number
    allele_number_numeric = min(allele_number_numerics)
    # we're going to rebuild the allele number from the components, with a leading zero if necessary
    lowest_allele = f"{allele_group}:{str(allele_number_numeric).zfill(2)}"

    return lowest_allele


def generate_allele_group_percentages(allele_group:Dict) -> Dict:
    """
    This function takes an allele group and returns a dictionary of the percentage of each allele in the allele group.

    Args:
        allele_group (Dict): An allele group dictionary.

    Returns:
        Dict: A dictionary of the percentage of each allele in the allele group.
    """
    # we'll initialize the percentages dictionary
    percentages = {}
    # we'll calculate the total number of pseudosequences in the allele group
    total = sum([allele_group['pseudosequences'][pseudosequence]['count'] for pseudosequence in allele_group['pseudosequences']])
    # now we'll iterate over the pseudosequences in the allele group
    for pseudosequence in allele_group['pseudosequences']:
        # we'll find the canonical allele for the pseudosequence (the one with the lowest allele number)
        canonical_allele = find_lowest_allele_number(allele_group['pseudosequences'][pseudosequence]['alleles'])
        # we'll calculate the percentage of the pseudosequence in the allele group
        percentages[slugify(canonical_allele)] = round(allele_group['pseudosequences'][pseudosequence]['count'] * 100/total, 3)
    # we'll sort the percentages dictionary by value in descending order
    percentages = dict(sorted(percentages.items(), key=lambda item: item[1], reverse=True))
    return percentages


def generate_allele_group_counts(allele_group:Dict) -> Dict:
    """
    This function takes an allele group and returns a dictionary of the count of each allele in the allele group.

    Args:
        allele_group (Dict): An allele group dictionary.
    
    Returns:
        Dict: A dictionary of the count of each allele in the allele group.
    """
    # we'll initialize the counts dictionary
    counts = {}
    # we'll iterate over the pseudosequences in the allele group
    for pseudosequence in allele_group['pseudosequences']:
        # we'll find the canonical allele for the pseudosequence (the one with the lowest allele number)
        canonical_allele = find_lowest_allele_number(allele_group['pseudosequences'][pseudosequence]['alleles'])
        # we'll add the count of the pseudosequence to the counts dictionary
        counts[slugify(canonical_allele)] = allele_group['pseudosequences'][pseudosequence]['count']
    # we'll sort the counts dictionary by value in descending order
    counts = dict(sorted(counts.items(), key=lambda item: item[1], reverse=True))
    return counts
    

def combine_allele_group_stats(counts:Dict, percentages:Dict) -> Dict:
    """
    This function takes a dictionary of counts and a dictionary of percentages and returns a combined dictionary.

    Note:
        The counts and percentages dictionaries must have the same keys/shape
    
    Args:
        counts (Dict): A dictionary of counts.
        percentages (Dict): A dictionary of percentages.

    Returns:
        Dict: A combined dictionary of counts and percentages.
    """
    # we'll initialize the combined dictionary
    combined = {}
    # we'll iterate over the counts
    for allele in counts:
        # we'll add the count and percentage to the combined dictionary
        combined[allele] = {'count':counts[allele], 'percentage':percentages[allele]}
    return combined


def top_n(group_stats:Dict, n:int=10) -> Tuple[List, List, List, List]:
    """
    This function takes a dictionary of counts and returns the top n labels and values.

    Args:
        counts (Dict): A dictionary of group stats including counts and percentages.
        n (int): The number of top values to return. Default is 10.

    Returns:
        Tuple[List, List, List, List]: A tuple containing the top n labels, top n percentages, top n counts, other labels, and other percentages.
    """
    # we'll sort the group stats by percentage in descending order
    sorted_data = sorted(group_stats.items(), key=lambda x: x[1]['percentage'], reverse=True)
    # we'll extract the top n labels and values
    i = 0
    top_n = []
    others = []
    for item in sorted_data:
        if i <= n:
            if item[1]['percentage'] > 5:
                top_n.append(item[0])
                i += 1
            else:
                others.append(item[0])
        else:
            others.append(item[0])
    labels = [item for item in top_n]
    percentages = [group_stats[item]['percentage'] for item in top_n]
    counts = [group_stats[item]['count'] for item in top_n]
    others_percentages = [group_stats[item]['percentage'] for item in others]
    return labels, percentages, counts, others, others_percentages


def deslugify_allele(text:str) -> str:
    """
    This function takes a slugified allele number and returns a deslugified allele number.

    TODO - this needs to be part of a centralised library of allele helper functions

    Args:
        text (str): A slugified allele number

    Returns:
        str: A deslugified allele number
    """
    components = text.split('_')
    allele_number = f"{components[0]}-{components[1]}*{components[2]}:{components[3]}".upper()
    return allele_number


def create_pie_chart(labels:List, percentages:List, counts:List, others:List, others_percentages:List, allele_group:str):
    """
    This function takes a list of labels and percentages and creates a pie chart.

    Args:
        labels (List): A list of labels.
        percentages (List): A list of percentages.
        counts (List): A list of counts.
        others (List): A list of other labels.
        others_percentages (List): A list of other percentages.
        allele_group (str): The allele group.
    
    Returns: 
        None
    """
    # first we need to rewrite the labels to include the percentage and to deslugify the allele number
    labels = [f"{deslugify_allele(label)} [{round(percentages[i], 1)}%]" for i, label in enumerate(labels)]

    # if there are other alleles, we'll add them to the labels and percentages
    if len(others) > 0:
        labels.append('Others')
        percentages.append(round(sum(others_percentages), 3))
    
    # we'll set the figure size
    figsize = 15
    
    # we'll create the figure object
    fig = Figure()
    # and set the size
    fig.set_figwidth(figsize+5)
    fig.set_figheight(figsize-3)
    ax = fig.subplots()
    bbox_props = dict(boxstyle="square,pad=0.2", fc="w", ec="k", lw=0)
    wedges, texts = ax.pie(percentages, textprops={'fontsize': 30}, wedgeprops=dict(width=0.3), startangle=-260, counterclock=False)

    kw = dict(arrowprops=dict(arrowstyle="-"), bbox=bbox_props, zorder=0, va="center")
    
    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = f"angle,angleA=0,angleB={ang}"
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(labels[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),horizontalalignment=horizontalalignment, **kw, size=30)

    # finally we'll save the figure as an SVG   
    filename = f"output/processed_data/pie_charts/allele_groups/{allele_group}.svg"
    fig.savefig(filename, format="svg", bbox_inches='tight', pad_inches=0.5)

    # and we'll close the figure object to free up memory
    fig = None
    pass


def create_allele_group_pie_chart(config:Dict, **kwargs):
    """
    This function takes a locus slug and creates a pie chart for each allele group.

    Args:
        locus_slug (str): The locus slug.

    Returns:
        None
    """
    # first, we'll load the pseudosequences for the locus
    locus = kwargs['locus']
    species_stem = kwargs['species_stem']
    
    locus_slug = f"{species_stem}_{locus.lower()}"

    input_filename = f"output/processed_data/pocket_pseudosequences/{locus_slug}.json"
    with open(input_filename, 'r') as f:
        pseudosequences = json.load(f)

    # next, we'll generate the allele groups from the pseudosequences
    allele_groups = generate_allele_groups(pseudosequences)

    # next, we'll iterate over the allele groups
    for allele_group in sorted(allele_groups.keys()):

        # we'll generate the allele group percentages and the counts
        allele_group_percentages = generate_allele_group_percentages(allele_groups[allele_group])
        allele_group_counts = generate_allele_group_counts(allele_groups[allele_group])

        # we'll combine the counts and percentages into a single dictionary
        allele_group_stats = combine_allele_group_stats(allele_group_counts, allele_group_percentages)

        # we'll get the top n labels and values
        labels, percentages, counts, others, others_percentages = top_n(allele_group_stats, 9)
        # and finally, we'll create the pie chart   
        create_pie_chart(labels, percentages, counts, others, others_percentages, allele_group)


def main():

    create_allele_group_pie_chart({}, locus='B', species_stem='hla')




if __name__ == '__main__':
    main()