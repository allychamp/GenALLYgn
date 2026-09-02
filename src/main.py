import glob
import os
import subprocess
import textwrap

import pandas as pd
import yaml

from src.blast_links import load_blast
from src.parser import parse_genbank
from src.plotter import plot_genomes
from src.distance import compute_genome_tree

with open("config.yaml") as f:
    config = yaml.safe_load(f)

# Paramètres
IDENTITY_THRESHOLD = config["blast"]["identity_threshold"]
COVERAGE_THRESHOLD = config["blast"]["coverage_threshold"]


# 1. Concatenate all protein fasta files
print("Step 1: Concatenating protein fasta files...")
faa_files = glob.glob(f"{config['paths']['faa_dir']}*.faa")
if not faa_files:
    faa_files = glob.glob(f"{config['paths']['faa_dir']}*.fasta")
with open(config["paths"]["all_proteins"], "w") as outfile:
    for faa in faa_files:
        with open(faa) as infile:
            outfile.write(infile.read())

# 2. Make BLAST database
print("Step 2: Building BLAST database...")
subprocess.run([
    "makeblastdb",
    "-in", config["paths"]["all_proteins"],
    "-dbtype", "prot"
], check=True)

# 3. Run all-vs-all BLAST
print("Step 3: Running all-vs-all BLAST...")
subprocess.run([
    "blastp",
    "-query", config["paths"]["all_proteins"],
    "-db", config["paths"]["all_proteins"],
    "-evalue", "1e-5",
    "-num_threads", str(config["blast"]["num_threads"]),
    "-outfmt", "6 qseqid sseqid pident length qlen slen bitscore evalue",
    "-out", config["paths"]["blast_file"]
], check=True)


# 4. Load gene colour lookup
cds_df = pd.read_csv(config["paths"]["cds_table"])
color_df = pd.read_csv(config["paths"]["color_table"])

# Strip Pharokka_ prefix to match the locus tags used everywhere else
cds_df["gene"] = cds_df["gene"].str.removeprefix("Pharokka_")

merged = cds_df.merge(color_df, left_on="category", right_on="function", how="left")
gene_color_map = dict(zip(merged["gene"], merged["color"]))


# 5. Parse the genomes
genomes = []
for gbk in sorted(glob.glob(f"{config['paths']['genomes_dir']}*.gbk")):
    parsed = parse_genbank(gbk)
    fname = os.path.basename(gbk).replace(".gbk", "")

    for i, g in enumerate(parsed):
        # One record per file is the normal case. If a file holds several, keep
        # the record name as a suffix so two genome objects never share a name:
        # the tree labels and the genome_y lookup in the plotter both key on it.
        g.name = fname if len(parsed) == 1 else f"{fname}_{g.name}"
        for gene in g.genes:
            gene.genome = g
        genomes.append(g)

print("Genomes parsed:", len(genomes))


# 6. Load BLAST hits
links = load_blast(
    config["paths"]["blast_file"],
    min_identity=IDENTITY_THRESHOLD,
    min_coverage=COVERAGE_THRESHOLD
)
print("Total links loaded:", len(links))
print("Thresholds — identity:", IDENTITY_THRESHOLD, "coverage:", COVERAGE_THRESHOLD)


# 7. Build the distance matrix and the tree.
# genome_names / distance_matrix / peq_matrix stay in the ORIGINAL genome order.
# Only `genomes` comes back permuted into plot order (bottom to top).
(genomes, linkage_matrix, distance_matrix, peq_matrix,
 genome_names) = compute_genome_tree(
    genomes, links,
    # Snap distances at or below this to exactly 0, so identical genomes join
    # as a flat line instead of a small bracket. Set to 0.0 to disable.
    zero_tolerance=config["plot"].get("zero_tolerance", 1e-3),
)


# 8. Map each gene to its genome
gene_index = {}
for genome in genomes:
    for gene in genome.genes:
        gene_index[gene.locus_tag] = (genome, gene)


# 9. Optional manual override of the row order.
# This breaks the correspondence between the tree and the rows, so drop the
# tree rather than draw one that contradicts the layout.
genome_order = config["plot"].get("genome_order", None)
if genome_order:
    name_to_genome = {g.name: g for g in genomes}
    missing = [n for n in genome_order if n not in name_to_genome]
    if missing:
        print("WARNING: genome_order names not found and skipped:", missing)
    dropped = [g.name for g in genomes if g.name not in set(genome_order)]
    if dropped:
        print("WARNING: genomes absent from genome_order and skipped:", dropped)

    genomes = [name_to_genome[name] for name in genome_order if name in name_to_genome]
    print("Manual genome_order in use — dendrogram disabled.")
    linkage_matrix = None


# 10. Best hit per gene per target genome
best_hits = {}
for hit in links:
    if hit.qseqid not in gene_index or hit.sseqid not in gene_index:
        continue
    g1, gene1 = gene_index[hit.qseqid]
    g2, gene2 = gene_index[hit.sseqid]
    if g1 is g2:
        continue
    key = (hit.qseqid, g2.name)
    if key not in best_hits or hit.identity > best_hits[key].identity:
        best_hits[key] = hit

# Convert to genome_links — the plotter filters these down to adjacent rows
genome_links = []
for hit in best_hits.values():
    g1, gene1 = gene_index[hit.qseqid]
    g2, gene2 = gene_index[hit.sseqid]
    genome_links.append((gene1, gene2, hit.identity))
print("Nombre de liens tracés :", len(genome_links))


# 11. Make every genome start at 0
for genome in genomes:
    genome.start = 0
    genome.end = genome.length


# 12. Plot
function_color_map = {
    textwrap.fill(func.capitalize(), width=config["plot"].get("legend_wrap", 25)): color
    for func, color in sorted(
        zip(color_df["function"], color_df["color"]),
        key=lambda x: x[0].lower()
    )
}

plot_cfg = config["plot"]

plot_genomes(
    genomes,
    genome_links,
    config["paths"]["output_svg"],
    identity_threshold=IDENTITY_THRESHOLD,
    gene_color_map=gene_color_map,
    function_color_map=function_color_map,
    spacing=plot_cfg["spacing"],
    linkage_matrix=linkage_matrix,
    # Original matrix order, NOT the reordered `genomes` list. scipy applies the
    # leaf permutation to `labels` itself.
    tree_labels=genome_names if linkage_matrix is not None else None,
    show_tree_labels=plot_cfg.get("show_tree_labels", True),
    show_genome_names=plot_cfg.get("show_genome_names", True),
    min_branch_display=plot_cfg.get("min_branch_display", 0.02),
    font_scale=plot_cfg.get("font_scale", 1.0),
    font_sizes=plot_cfg.get("fonts", None),
    tree_gap=plot_cfg.get("tree_gap", None),
)
print("Wrote", config["paths"]["output_svg"])