# src/plotter.py

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram
import numpy as np


# Defaults, overridable per key from config.yaml (plot.fonts).
DEFAULT_FONT_SIZES = {
    "genome_name": 18,    # names beside each genome line
    "tree_label": 14,     # names at the branch tips
    "tree_title": 16,     # "Tree"
    "tree_axis": 12,      # distance scale under the tree
    "position_axis": 16,  # "Position (bp)" and its ticks
    "colorbar": 15,       # "Identity (%)" and its ticks
    "legend": 15,         # gene function legend entries
    "legend_title": 17,   # "Gene function"
}


def draw_gene(ax, gene, y, height=0.4, color="#5a5a5f"):
    if gene.strand == 1:
        x1 = gene.start
        x2 = gene.end
    else:
        x1 = gene.end
        x2 = gene.start

    length = abs(x2 - x1)
    head = min(length * 0.8, 1000)

    if gene.strand == 1:
        points = [
            (x1, y - height / 3),
            (x2 - head, y - height / 3),
            (x2, y),
            (x2 - head, y + height / 3),
            (x1, y + height / 3),
        ]
    else:
        points = [
            (x1, y - height / 3),
            (x2 + head, y - height / 3),
            (x2, y),
            (x2 + head, y + height / 3),
            (x1, y + height / 3),
        ]

    ax.add_patch(
        Polygon(points, closed=True, facecolor=color, edgecolor="black", linewidth=1)
    )


def draw_link(ax, gene1, gene2, y1, y2, color, alpha):
    poly = [
        (gene1.start, y1),
        (gene1.end, y1),
        (gene2.end, y2),
        (gene2.start, y2),
    ]
    ax.add_patch(
        Polygon(poly, closed=True, facecolor=color, edgecolor=None, alpha=alpha, zorder=0)
    )


def plot_genomes(genomes, links, output_file, identity_threshold=10,
                 gene_color_map=None, function_color_map=None,
                 spacing=1.0, linkage_matrix=None, tree_labels=None,
                 show_tree_labels=True, show_genome_names=True,
                 min_branch_display=0.02, font_scale=1.0, font_sizes=None,
                 tree_gap=None):
    """
    `genomes` must already be in plot order (bottom to top).

    `tree_labels` must be the genome names in the ORIGINAL distance-matrix
    order, NOT the plot order. scipy's dendrogram applies the leaf permutation
    to `labels` itself; handing it pre-permuted names applies the permutation
    twice and the tree stops lining up with the rows.

    `show_tree_labels` prints the names at the branch tips. Useful for checking
    that the tree and the rows agree; set it False once you trust the output.
    `show_genome_names` draws the names beside each genome line. Turning one of
    the two off avoids printing every name twice.

    `font_scale` multiplies every text size. `font_sizes` overrides individual
    ones by name; see DEFAULT_FONT_SIZES for the keys. Both come from
    config.yaml via main.py.
    """

    sizes = dict(DEFAULT_FONT_SIZES)
    if font_sizes:
        unknown = set(font_sizes) - set(DEFAULT_FONT_SIZES)
        if unknown:
            raise ValueError(
                f"Unknown font size key(s) in config: {sorted(unknown)}. "
                f"Valid keys: {sorted(DEFAULT_FONT_SIZES)}"
            )
        sizes.update(font_sizes)

    FS_GENOME_NAME = sizes["genome_name"] * font_scale
    FS_TREE_LABEL = sizes["tree_label"] * font_scale
    FS_TREE_TITLE = sizes["tree_title"] * font_scale
    FS_TREE_AXIS = sizes["tree_axis"] * font_scale
    FS_POSITION = sizes["position_axis"] * font_scale
    FS_CBAR = sizes["colorbar"] * font_scale
    FS_LEGEND = sizes["legend"] * font_scale
    FS_LEGEND_TITLE = sizes["legend_title"] * font_scale

    fig_width = 20
    fig = plt.figure(figsize=(fig_width, 4 + len(genomes) * spacing))

    if linkage_matrix is not None:
        # The gap between the panels has to hold two sets of labels: the tree
        # leaf names, which overhang ax_tree to the right, and the genome names,
        # which are drawn left of x=0 and so overhang the alignment axes to the
        # left. Both grow with font_scale, so a fixed wspace lets them collide.
        if tree_gap is None:
            # ~0.6 em per character is a decent estimate for the default face.
            char_w = 0.6 / 72.0
            label_in = 0.0
            if show_tree_labels and tree_labels:
                label_in += max(len(str(t)) for t in tree_labels) * char_w * FS_TREE_LABEL
            if show_genome_names:
                label_in += max(len(g.name) for g in genomes) * char_w * FS_GENOME_NAME
            # Padding scales too: the per-character estimate drifts low as the
            # text grows, so a fixed pad runs out at large font_scale.
            gap_in = label_in + 0.3 * max(font_scale, 1.0)

            # wspace is a fraction of the mean axes width, and the axes share
            # whatever the gap leaves behind.
            usable_in = fig_width * (1.0 - 0.025)
            wspace = 2 * gap_in / max(usable_in - gap_in, 1.0)
            wspace = float(np.clip(wspace, 0.15, 3.0))
        else:
            wspace = tree_gap

        gs = fig.add_gridspec(1, 2, width_ratios=[0.75, 5], wspace=wspace)
        ax_tree = fig.add_subplot(gs[0, 0])
        ax = fig.add_subplot(gs[0, 1])
    else:
        ax = fig.add_subplot(111)
        ax_tree = None

    ordered_genomes = genomes

    # ======================
    # 0) Dendrogram
    # ======================
    if linkage_matrix is not None:
        if tree_labels is None:
            raise ValueError(
                "tree_labels is required when linkage_matrix is given. Pass the "
                "genome names in original matrix order (the fifth return value "
                "of compute_genome_tree)."
            )
        if len(tree_labels) != len(ordered_genomes):
            raise ValueError(
                f"tree_labels has {len(tree_labels)} entries but there are "
                f"{len(ordered_genomes)} genomes."
            )

        # A merge at distance 0 has no arms: scipy draws each link as
        # (0, y1) -> (h, y1) -> (h, y2) -> (0, y2), so h = 0 leaves only the
        # vertical crossbar sitting on the leaf tips. Give the drawing a floor
        # so those links render as a normal bracket. This copies the linkage
        # matrix, so distances, PEQ and the leaf order are untouched -- only
        # the rendered geometry changes. Set min_branch_display=0.0 to see the
        # true zero-length branches.
        Z_draw = np.asarray(linkage_matrix, dtype=float).copy()
        if min_branch_display > 0 and Z_draw[:, 2].max() > 0:
            floor = min_branch_display * Z_draw[:, 2].max()
            Z_draw[:, 2] = np.maximum(Z_draw[:, 2], floor)

        d = dendrogram(
            Z_draw,
            orientation="left",
            labels=list(tree_labels),
            ax=ax_tree,
            color_threshold=0,
            above_threshold_color="black",
        )

        # Sanity check: the leaves, bottom to top, must be the plotted rows,
        # bottom to top. If this fires, the genome order and the tree disagree
        # (e.g. config genome_order overrode the clustering).
        if d["ivl"] != [g.name for g in ordered_genomes]:
            print(
                "WARNING: dendrogram leaf order does not match genome plot order.\n"
                f"  tree: {d['ivl']}\n"
                f"  rows: {[g.name for g in ordered_genomes]}"
            )

        ax_tree.set_title("Proteomic Equivalence \n Quotient Tree", fontsize=FS_TREE_TITLE)
        ax_tree.spines["top"].set_visible(False)
        ax_tree.spines["right"].set_visible(False)
        ax_tree.spines["left"].set_visible(False)

        # Distance scale. With orientation="left" the root sits on the left, so
        # the axis runs from max distance on the left down to 0 at the leaves.
        ax_tree.spines["bottom"].set_visible(True)
        ax_tree.tick_params(axis="x", labelsize=FS_TREE_AXIS)
        ax_tree.set_xlabel("Distance (1 - PEQ)", fontsize=FS_TREE_AXIS)

        if show_tree_labels:
            # Put the names at the branch tips (right edge), so each label sits
            # on the same row as the genome it belongs to.
            ax_tree.yaxis.tick_right()
            ax_tree.tick_params(axis="y", length=0, labelsize=FS_TREE_LABEL)
        else:
            ax_tree.set_yticks([])

    # ======================
    # Rebuild adjacent pairs based on genome order
    # ======================
    adjacent_pairs = set()
    for i in range(len(ordered_genomes) - 1):
        pair = tuple(sorted([ordered_genomes[i].name, ordered_genomes[i + 1].name]))
        adjacent_pairs.add(pair)

    links = [
        (gene1, gene2, identity)
        for gene1, gene2, identity in links
        if tuple(sorted([gene1.genome.name, gene2.genome.name])) in adjacent_pairs
    ]

    # ======================
    # Colormap setup
    # ======================
    cmap = LinearSegmentedColormap.from_list("white_to_black", ["white", "black"])
    norm = mcolors.Normalize(vmin=identity_threshold, vmax=100)

    # ======================
    # 1) Draw genome lines
    # ======================
    genome_y = {}
    y = 0
    for genome in ordered_genomes:
        genome_y[genome.name] = y
        ax.hlines(y, genome.start, genome.end, linewidth=4, color="gray", zorder=0)
        if show_genome_names:
            ax.text(genome.start - 0.02 * genome.end, y, genome.name,
                    ha="right", va="center", fontsize=FS_GENOME_NAME, weight="bold")
        y += spacing

    # ======================
    # 2) Draw genes
    # ======================
    for genome in ordered_genomes:
        y = genome_y[genome.name]
        for gene in genome.genes:
            color = gene_color_map.get(gene.locus_tag, "#5a5a5f") if gene_color_map else "#5a5a5f"
            draw_gene(ax, gene, y, color=color)

    # ======================
    # 3) Draw links
    # ======================
    for gene1, gene2, identity in links:
        y1 = genome_y[gene1.genome.name]
        y2 = genome_y[gene2.genome.name]
        if identity > identity_threshold:
            color = mcolors.to_hex(cmap(norm(identity)))
        else:
            color = "none"
        draw_link(ax, gene1, gene2, y1, y2, color, alpha=1)

    # ======================
    # 4) Formatting
    # ======================
    n = len(ordered_genomes)
    ax.set_ylim(-spacing * 0.5, (n - 1) * spacing + spacing * 0.5)
    ax.set_xlabel("Position (bp)", fontsize=FS_POSITION)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    # scipy puts leaf k at y = 10k + 5 with ylim (0, 10n), so row k sits at the
    # same fractional height as leaf k. Set it explicitly rather than relying on
    # the default, so the panels stay locked together.
    if ax_tree is not None:
        ax_tree.set_ylim(0, 10 * n)

    # ======================
    # 5) Colorbar
    # ======================
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    plt.subplots_adjust(left=0.025, right=1, top=0.95, bottom=0.12)
    fig.canvas.draw()

    x_start_display = ax.transData.transform((min(g.start for g in ordered_genomes), 0))[0]
    fig_width_display = fig.get_figwidth() * fig.dpi
    x0_fig = x_start_display / fig_width_display

    cax = fig.add_axes([x0_fig, 0.02, 0.2, 0.02])
    cbar = plt.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.set_label("Identity (%)", fontsize=FS_CBAR)
    cbar.ax.tick_params(labelsize=FS_CBAR)
    ax.tick_params(axis="x", labelsize=FS_POSITION)
    cbar.set_ticks(np.arange(identity_threshold, 100.01, 10))

    # ======================
    # 6) Gene function legend
    # ======================
    if function_color_map:
        handles = [Patch(facecolor=color, edgecolor="black", linewidth=0.6, label=func)
                   for func, color in function_color_map.items()]
        ax.legend(handles=handles, title="Gene function",
                  loc="upper left",
                  bbox_to_anchor=(1.01, 1),
                  bbox_transform=ax.transAxes,
                  frameon=True, fontsize=FS_LEGEND, title_fontsize=FS_LEGEND_TITLE)

    plt.subplots_adjust(left=0.025, right=1, top=0.95, bottom=0.1)
    plt.savefig(output_file, dpi=600, bbox_inches="tight")
    plt.close()