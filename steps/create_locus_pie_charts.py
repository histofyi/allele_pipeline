from typing import Dict, List, Tuple

import base64
from io import BytesIO

from matplotlib.figure import Figure
import numpy as np

import json

def top_n(dataset:Dict, n:int=10):
    # Sort the dictionary by percentage in descending order
    sorted_data = sorted(dataset.items(), key=lambda x: x[1]['percent'], reverse=True)
    
    # Extract the top 10 labels and values
    labels = [item[0] for item in sorted_data[:n]]
    values = [item[1]['percent'] for item in sorted_data[:n]]
    
    # Extract the rest of the labels and values for others category
    others = [item[0] for item in sorted_data[n:]]
    others_values = [item[1]['percent'] for item in sorted_data[n:]]
    
    return labels, values, others, others_values

def deslugify_allele_group(text:str) -> str:
    elements = text.split('_')
    return f"{elements[0]}-{elements[1]}*{elements[2]}".upper()


def generate_allele_group_pie_chart(allele_groups:Dict, allele_count:int, locus:str) -> Tuple[str, str, str]:
    labels = []
    values = []
    others = []

    for allele_group in allele_groups:
        allele_groups[allele_group]['percent'] = allele_groups[allele_group]['count']*100/allele_count

    labels, values, others, others_values = top_n(allele_groups, 9)

    labels = [f"{deslugify_allele_group(allele_group)} [{round(allele_groups[allele_group]['percent'], 1)}%]" for allele_group in labels]

    others_percent = sum(others_values)

    if len(others) > 0:
        labels.append('Others')
        values.append(others_percent)
    figsize = 15

    fig = Figure()
    fig.set_figwidth(figsize+4)
    fig.set_figheight(figsize-3)
    ax = fig.subplots()
    bbox_props = dict(boxstyle="square,pad=0.2", fc="w", ec="k", lw=0)
    wedges, texts = ax.pie(values, textprops={'fontsize': 30}, wedgeprops=dict(width=0.3), startangle=-260, counterclock=False)

    kw = dict(arrowprops=dict(arrowstyle="-"), bbox=bbox_props, zorder=0, va="center")
    
    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = f"angle,angleA=0,angleB={ang}"
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(labels[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),horizontalalignment=horizontalalignment, **kw, size=32)

    # Save it to a temporary buffer.
    png = BytesIO()
    svg = BytesIO()


    filestem = f"output/processed_data/pie_charts/{locus.lower()}"
    
    fig.savefig(png, format="png")
    fig.savefig(f"{filestem}.png", format="png")
    
    fig.savefig(svg, format="svg")
    fig.savefig(f"{filestem}.svg", format="svg")

    png_data = base64.b64encode(png.getbuffer()).decode("ascii")    
    svg_data = base64.b64encode(svg.getvalue()).decode("ascii")

    with open(f"{filestem}_svg.txt", 'w') as svg_file:
        svg_file.write(svg_data)
    with open(f"{filestem}_png.txt", 'w') as png_file:
        png_file.write(png_data)

    return png_data, svg_data, None



def create_locus_pie_charts(config:Dict, **kwargs) -> None:
    """
    This function takes allele group lists for a locus and makes a pie chart of the top allele groups.

    Args:
        locus (str): the locus to be parsed
        verbose (bool): whether specific information is output to the terminal, for large sequence sets this can be overwhelming and significantly slow down the function
    Returns:

    """
    locus = kwargs['locus']
    species_stem = kwargs['species_stem']
    
    locus_slug = f"{species_stem}_{locus.lower()}"


    input_filename = f"output/processed_data/allele_groups/{locus_slug}.json"

    with open(input_filename, 'r') as allele_groups_file:
        allele_groups = json.load(allele_groups_file)
    

    allele_count = 0
    allele_group_count = 0

    allele_group_stats = {}

    for allele_group in allele_groups:
        alleles_count = len(allele_groups[allele_group])
        allele_group_stats[allele_group] = {'count':alleles_count, 'percent':0}

        allele_count += alleles_count
        allele_group_count += 1


    png_data, svg_data, alt_text = generate_allele_group_pie_chart(allele_group_stats, allele_count, locus)

    pass
