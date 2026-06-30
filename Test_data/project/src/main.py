import pandas as pd
from src.blast_links import load_blast
from src.parser import parse_genbank
from src.plotter import plot_genomes
from collections import defaultdict
import textwrap
import glob
import os
# ⚙️ Parameters
IDENTITY_THRESHOLD = 85



# Load gene color lookup
cds_df = pd.read_csv("/home/champa/DATA/PHAGE_genome_analysis/Champoux_A/Analysis_05-2026/Analysis/clinker_alignements/merged_cds_predictions.csv")
color_df = pd.read_csv("/home/champa/DATA/PHAGE_genome_analysis/Champoux_A/Analysis_05-2026/Analysis/clinker_alignements/function_color_mapping.csv")

# Strip Pharokka_ prefix to match your locus tags
cds_df["gene"] = cds_df["gene"].str.removeprefix("Pharokka_")

# Merge to get gene -> color
merged = cds_df.merge(color_df, left_on="category", right_on="function", how="left")
gene_color_map = dict(zip(merged["gene"], merged["color"]))

def order_genomes_by_similarity(genomes, genome_links):
    """
    A function to class genomes by similarity to display the closest together in the graph 
    genomes: a list of genome objects
    genome_links: a list of (gene1, gene2, identity) tuples representing connections between genes from different genomes, with a similarity score
    """
    
    #Defying two dictionnaries to store the sum of each similarity score for each genome pair (pair_scores) and 
    #how many links exist between each pair (pair_counts)
    pair_scores = defaultdict(float)
    pair_counts = defaultdict(int)

    #Iterating over every connections between genomes
    for gene1, gene2, identity in genome_links:
        #Creating pairs of each genome combination
        key = tuple(sorted([gene1.genome.name, gene2.genome.name]))
        # Sum gene identity into the pair_score value associated with the key (genome pair)
        pair_scores[key] += identity
        # Add one the the count of gene links (value) between the key (genome pair)
        pair_counts[key] += 1
    #Then the pair score becomes the average of similarity of all genes in the key (genome pair)
    for key in pair_scores:
        pair_scores[key] /= pair_counts[key]

    if not pair_scores:
        print("WARNING: no pairs found, keeping original order")
        return genomes
    #Extracting in genome_names variable just the name strings from the genome objects
    genome_names = [g.name for g in genomes]
    #Finding the genome pair with the highest average similarity and uses them as the starting two elements of the chain.
    #This chain will contain all the genomes names in the right order. The genomes not in the ''best_pair'' will be store
    # in the remaining variable. 
    best_pair = max(pair_scores, key=pair_scores.get)
    chain = list(best_pair)
    remaining = [n for n in genome_names if n not in chain]

    #The loop will continue until until all remainning has been placed in the chain
    while remaining:
        best_score = -1
        best_genome = None
        best_end = None
    #Checking how similar the next genome in remainning is to the left end of the chain. 
    # Uses 0 if no link exists between them.
        for name in remaining:
            left_key = tuple(sorted([chain[0], name]))
            left_score = pair_scores.get(left_key, 0)
            #Same for the right end of the chain
            right_key = tuple(sorted([chain[-1], name]))
            right_score = pair_scores.get(right_key, 0)
        #If the similarity is higher than previous matches, the new match will replace the best candidates will be updated
            if left_score > best_score:
                best_score = left_score
                best_genome = name
                best_end = "left"
            if right_score > best_score:
                best_score = right_score
                best_genome = name
                best_end = "right"

        if best_end == "left":
            chain.insert(0, best_genome)
        else:
            chain.append(best_genome)
        remaining.remove(best_genome)
#Place the name back into a dictionnary
    name_to_genome = {g.name: g for g in genomes}
    return [name_to_genome[name] for name in chain] 

#Use the function create in the parser.py to return genomes object 
genomes = []
# Iterates over all the gbk in genomes file
for gbk in glob.glob("data/genomes/*.gbk"):
    parsed = parse_genbank(gbk)
    #Attributing the file name as genome's name
    fname = os.path.basename(gbk).replace('.gbk', '')

    # For every gene inside this genome, setting the genome attributes of the gene object (define in parser)
    for g in parsed:
        g.name = fname
        for gene in g.genes:
            gene.genome = g 
        genomes.append(g)  


# gene_lookup = {}

# for genome in genomes:
#     for gene in genome.genes:
#         gene_lookup[gene.locus_tag] = gene

# print("Nombre de gènes :", len(gene_lookup))

#Use function define in blast_link
links = load_blast("data/blast/allPV_vs_allPV.tsv", min_identity=IDENTITY_THRESHOLD)

# for tag in list(gene_lookup.keys())[:5]:
#     print(f"  {tag}")
# print("Nombre de hits BLAST :", len(links))
# print("Sample locus tags:", list(gene_lookup.keys())[:5])
# print("Sample BLAST IDs:", [links[i].qseqid for i in range(5)])

#Creating a dictionnary (gene_index) to map each gene with the appropriate genomes (could probably have been done with the attributes, tcheck that later)
gene_index = {}
for genome in genomes:
    for gene in genome.genes:
        gene_index[gene.locus_tag] = (genome, gene)

# 
all_genome_links = []
for hit in links:
    #Make sure to only keep hits in gbk and blast 
    if hit.qseqid not in gene_index or hit.sseqid not in gene_index:
        continue
    #g1 is defined as genome of the query gene (gene1), g2 is defined as genome of the subject gene (g2)
    g1, gene1 = gene_index[hit.qseqid]
    g2, gene2 = gene_index[hit.sseqid]
    #Remove self hit (as precaution because thechnically is as already be removed)
    if g1 == g2:
        continue
    all_genome_links.append((gene1, gene2, hit.identity))

#Using function defined earlier to obtain the order genomes should appear in plot
genomes = order_genomes_by_similarity(genomes, all_genome_links)

# Making set of adjacent genome only
adjacent_pairs = set()
for i in range(len(genomes) - 1):
    pair = tuple(sorted([genomes[i].name, genomes[i+1].name]))
    adjacent_pairs.add(pair)

#Finding best hit in adjacent genomes only. Best hits will be stored in a dictionnary
best_hits = {}
#Iterating trought blast results (again)
for hit in links:
    if hit.qseqid not in gene_index or hit.sseqid not in gene_index:
        continue
    g1, gene1 = gene_index[hit.qseqid]
    g2, gene2 = gene_index[hit.sseqid]
    if g1 == g2:
        continue
    pair = tuple(sorted([g1.name, g2.name]))
    #Keeping only the hit that are between genomes in adjacent_pairs set
    if pair not in adjacent_pairs:
        continue
    # Key is (query_gene, adjacent_genome) to get best hit per neighbor
    key = (hit.qseqid, g2.name)
    #Updating key to only keep the best hit
    if key not in best_hits or hit.identity > best_hits[key].identity:
        best_hits[key] = hit

# Converting the best hit into visual links on the figure
genome_links = []
for hit in best_hits.values():
    g1, gene1 = gene_index[hit.qseqid]
    g2, gene2 = gene_index[hit.sseqid]
    genome_links.append((gene1, gene2, hit.identity))

print("Nombre de liens tracés :", len(genome_links))

# Assuring that all genomes start at 0
for genome in genomes:
    genome.start = 0
    genome.end = genome.length

# Plotting the figure with the function define in the plotter.py

function_color_map = {
    textwrap.fill(func.capitalize(), width=25): color
    for func, color in sorted(
        zip(color_df["function"], color_df["color"]),
        key=lambda x: x[0].lower()
    )
}
plot_genomes(genomes, genome_links, "genome_comparison_all_PV.svg",
             identity_threshold=IDENTITY_THRESHOLD,
             gene_color_map=gene_color_map,
             function_color_map=function_color_map,
             spacing=0.8)