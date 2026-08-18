# src/plotter.py

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory

  
import numpy as np


#Definition of a function that will represent each gene in an arrow. The gene has been define in the parser.py, 
# The height is the tickness of the arrow of eahc gene (define by user). The color will change
def draw_gene(ax, gene, y, height=0.4, color="#5a5a5f"):
    # function_to_color = {
    # 'Tail': '#dcf3ff',
    # 'Collars': '#baf2ef',
    # 'Head-tail joinning':'#0013de',
    # 'Portal': '#396d7c',
    # 'Capsid': '#4e4cb0',

    # 'Helicase': "#ffaa1d",
    # 'Nuclease': "#ffdf00",
    # 'Terminase': "#ffff00",
    # 'Integration': "#b18914",
    # 'DNA polymerase': '#ffe05f',
    # 'Annealing': "#f7a151",
    # 'Primase': "#fff380",
    # 'Replication initiation': "#e77500",

    # 'Activator': '#ba0012',
    # 'Repressor': "#a04848",

    # 'Toxin/anti-toxin': "#c2f74f",
    # 'Anti-restriction': '#74d600',
    # 'Sir2': '#028900',
    # 'CRISPR/anti-CRISPR': '#13762e',
    # 'Super-infection exclusion': '#cce451',

    # 'Endolysin': '#fa638e',
    # 'Spanin': '#ffabc8',
    # 'Lysis inhibition': '#fcdfe3',
    # 'Holin': '#ff3377',

    # 'Packaging/assembly': '#bfbaf2',
    # 'Ejection': "#4a1170",
    # 'Reductase': "#afa9a9",
    # 'Phosphorylation': "#20706c",
    # 'RNA-associated': "#681b1b",
    # 'Nucleotide metabolism': "#5e4120",
    # 'Transferase': "#9c755b",
    # 'Cell wall depolymerase': "#c9c6ec",
    # 'Adsorption-related': "#b48ca8",
    # 'Unknown': "#5a5a5f"
    # }
    #That flip the gene if it's identify on the antisens strand
    if gene.strand == 1:
        x1 = gene.start
        x2 = gene.end
    else:
        x1 = gene.end
        x2 = gene.start
    # function_color = ()
    # for f,c in function_color.items():
    #     if f in gene.function:
    #         function_color = function_to_color.keys()
    

    # Lenght is define as the lenght of the gene and head is the head of the arrow that will be use later
    length = abs(x2 - x1)
    head = min(length*0.8, 1000)

    # Then points are defined to draw the arrow later. 
    if gene.strand == 1:
        points = [
            (x1, y - height/3),
            (x2 - head, y - height/3),
            (x2, y),
            (x2 - head, y + height/3),
            (x1, y + height/3)
        ]
    else:
        points = [
            (x1, y - height/3),
            (x2 + head, y - height/3),
            (x2, y),
            (x2 + head, y + height/3),
            (x1, y + height/3)
        ]

    # The arrow are added 
    ax.add_patch(
        Polygon(
            points,
            closed=True,
            facecolor=color,  # ✅ use passed color
            edgecolor="black",
            linewidth=1
        )
    )

# Here is the function to draw the link between the genomes
def draw_link(ax, gene1, gene2, y1, y2, color, alpha):
    poly = [
        (gene1.start, y1),
        (gene1.end,   y1),
        (gene2.end,   y2),
        (gene2.start, y2)
    ]
   
    ax.add_patch(
        Polygon(poly, closed=True, facecolor=color, edgecolor=None, alpha=alpha, zorder =0)
    )

# Then the function to draw the plot
def plot_genomes(genomes, links, output_file, identity_threshold=10, 
                 gene_color_map=None, function_color_map=None,spacing=1.0):
    fig, ax = plt.subplots(figsize=(20,4+len(genomes)*spacing))

    # Colormap setup
    
    cmap = LinearSegmentedColormap.from_list("white_to_black", ["white", "black"])
    norm = mcolors.Normalize(vmin=identity_threshold, vmax=100)

    # 1) Draw the lines representing each genomes
    # genome_y wil contain all the y coordinates for each genomes 
    genome_y = {}
    y = 0
    # Then drawing the gray line that will be behind each arrow
    for genome in genomes:
        genome_y[genome.name] = y
        ax.hlines(y, genome.start, genome.end, linewidth=4, color="gray", zorder=0)
        # 
        ax.text(genome.start - 0.02 * genome.end, y, genome.name,
                ha="right", va="center",fontsize=14, weight = 'bold' )
        y += spacing

    # 2)
    for genome in genomes:
        y = genome_y[genome.name]
        for gene in genome.genes:
            color = gene_color_map.get(gene.locus_tag, "#5a5a5f") if gene_color_map else "#5a5a5f"
            draw_gene(ax, gene, y, color=color)

    # ======================
    # 3) Dessin des liens
    # ======================
    # print(links)
    for gene1, gene2, identity in links:
        y1 = genome_y[gene1.genome.name] 
        y2 = genome_y[gene2.genome.name]
        if identity > identity_threshold:
            color = mcolors.to_hex(cmap(norm(identity)))
        else:
            color = "none"
        alpha = 1

        draw_link(ax, gene1, gene2 , y1, y2, color, alpha)

    # ======================
    # 4) Mise en forme
    # ======================
    ax.set_ylim(-spacing * 0.5, (len(genomes) - 1) * spacing + spacing * 0.5)
    ax.set_xlabel("Position (bp)", fontsize = 12)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
# 5) Colorbar
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    plt.subplots_adjust(left=0.025, right=1, top=0.95, bottom=0.12)

    fig.canvas.draw() 
    x_start_display = ax.transData.transform((min(g.start for g in genomes), 0))[0]
    fig_width_display = fig.get_figwidth() * fig.dpi
    x0_fig = x_start_display / fig_width_display

    cax = fig.add_axes([x0_fig, 0.02, 0.2, 0.02])
    cbar = plt.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.set_label("Identity (%)", fontsize=12)
    cbar.ax.tick_params(labelsize=12)
    ax.tick_params(axis="x", labelsize=11)
  
    cbar.set_ticks(np.arange(identity_threshold, 100.01, 10)) 
    # 6) Gene function legend
    if function_color_map:
        handles = [Patch(facecolor=color, edgecolor="black", linewidth=0.6, label=func)
                for func, color in function_color_map.items()]

        ax.legend(handles=handles, title="Gene function",
                loc="upper left",
                bbox_to_anchor=(1.01, 0),        # ← right side, bottom
                bbox_transform=ax.transAxes,
                frameon=True, fontsize=12, title_fontsize=14)
    plt.subplots_adjust(left=0.025, right=1, top=0.95, bottom=0.1)
    plt.savefig(output_file, dpi=600, bbox_inches="tight")
    plt.close()