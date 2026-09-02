# src/distance.py

from collections import defaultdict

import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram, optimal_leaf_ordering
from scipy.spatial.distance import squareform

# BLAST protein alignment lengths are in amino acids; gene.length is in bp.
AA_TO_BP = 3


def build_gene_index(genomes):
    """
    locus_tag -> (genome, gene)

    Warns on duplicate locus tags, since a collision silently reassigns a gene
    to the wrong genome and corrupts every pair it takes part in.
    """
    gene_index = {}
    duplicates = set()
    for genome in genomes:
        for gene in genome.genes:
            if gene.locus_tag in gene_index:
                duplicates.add(gene.locus_tag)
            gene_index[gene.locus_tag] = (genome, gene)

    if duplicates:
        print(
            f"WARNING: {len(duplicates)} duplicate locus tag(s) across genomes, "
            f"e.g. {sorted(duplicates)[:5]}. Hits for these are unreliable."
        )
    return gene_index


def compute_peq_matrix(genomes, hits):
    """
    PEQ = AF * AAI for every genome pair.

    AF  = (aligned bp in A + aligned bp in B) / (total CDS bp in A + total CDS bp in B)
    AAI = alignment-length-weighted mean identity over the same best hits

    Best hits are taken per direction: for each gene of A its best hit in B, and
    independently for each gene of B its best hit in A. This is what stops a
    single gene from being counted once per gene that hits it.

    Returns (genome_names, peq) where genome_names is in the *input* order of
    `genomes`, which is also the row/column order of the matrix.
    """
    gene_index = build_gene_index(genomes)

    genome_names = [g.name for g in genomes]
    # Key on object identity, not name: two records in one .gbk can share a name.
    genome_to_idx = {id(g): i for i, g in enumerate(genomes)}

    n = len(genomes)
    peq = np.eye(n)

    if len(set(genome_names)) != n:
        print(
            "WARNING: duplicate genome names. The matrix is still correct, but "
            "labels in the tree will be ambiguous."
        )

    total_length = [sum(gene.length for gene in g.genes) for g in genomes]

    # best[(i, j)][locus_tag_in_i] = (identity, aligned_bp)
    best = defaultdict(dict)

    for hit in hits:
        q = gene_index.get(hit.qseqid)
        s = gene_index.get(hit.sseqid)
        if q is None or s is None:
            continue

        genome_q, gene_q = q
        genome_s, gene_s = s

        i = genome_to_idx[id(genome_q)]
        j = genome_to_idx[id(genome_s)]
        if i == j:
            continue

        aligned_bp = hit.length * AA_TO_BP

        # Register both directions from the single hit. All-vs-all BLAST usually
        # reports both, but load_blast only filters on query coverage, so one
        # direction can be dropped. Doing it here keeps the matrix symmetric.
        for (a, b, locus, gene) in (
            (i, j, hit.qseqid, gene_q),
            (j, i, hit.sseqid, gene_s),
        ):
            # A local alignment cannot cover more of the gene than the gene has.
            covered = min(aligned_bp, gene.length)
            prev = best[(a, b)].get(locus)
            if prev is None or hit.identity > prev[0]:
                best[(a, b)][locus] = (hit.identity, covered)

    for i in range(n):
        for j in range(i + 1, n):
            forward = best.get((i, j))
            reverse = best.get((j, i))
            if not forward and not reverse:
                continue

            shared_i = sum(v[1] for v in forward.values()) if forward else 0.0
            shared_j = sum(v[1] for v in reverse.values()) if reverse else 0.0

            denom = total_length[i] + total_length[j]
            if denom == 0:
                continue

            AF = min((shared_i + shared_j) / denom, 1.0)

            weight = shared_i + shared_j
            if weight == 0:
                continue

            weighted_identity = 0.0
            for side in (forward, reverse):
                if side:
                    weighted_identity += sum(
                        (v[0] / 100.0) * v[1] for v in side.values()
                    )
            AAI = weighted_identity / weight

            score = AF * AAI
            peq[i, j] = score
            peq[j, i] = score

    return genome_names, peq


def compute_distance_matrix(peq_matrix, zero_tolerance=0.00000):
    """
    zero_tolerance: distances at or below this value are snapped to exactly 0.

    Two genomes only reach distance 0 when AF and AAI are both exactly 1.0. One
    CDS at 99.8% identity, or one alignment that stops short of the stop codon,
    leaves a residue of ~1e-4. In the dendrogram that residue is drawn as a
    short bracket at the leaf tips instead of a flat line, because scipy draws
    each link as (0, y1) -> (h, y1) -> (h, y2) -> (0, y2) and only collapses to
    a straight line when h is exactly 0.

    This changes the reported distances, not just the drawing, so pick the value
    deliberately: 1e-3 means "PEQ above 0.999 counts as identical".
    """
    dist = 1.0 - peq_matrix
    # Guard against tiny negatives from PEQ slightly exceeding 1.0.
    np.clip(dist, 0.0, None, out=dist)
    # Force exact symmetry so squareform is happy.
    dist = (dist + dist.T) / 2.0

    if zero_tolerance > 0:
        dist[dist <= zero_tolerance] = 0.0

    np.fill_diagonal(dist, 0.0)
    return dist


def build_tree(distance_matrix, method="average", optimal_order=True):
    """
    Returns (Z, order) where `order` is the left-to-right / bottom-to-top leaf
    order, expressed as indices into the *original* matrix order.
    """
    n = distance_matrix.shape[0]
    if n < 2:
        return None, list(range(n))

    condensed = squareform(distance_matrix, checks=False)
    Z = linkage(condensed, method=method)

    # This is what actually puts similar genomes next to each other. Plain
    # linkage fixes the topology but leaves each node free to flip either way.
    if optimal_order and n > 2:
        Z = optimal_leaf_ordering(Z, condensed)

    order = dendrogram(Z, no_plot=True)["leaves"]
    return Z, order


def order_genomes(genomes, order):
    return [genomes[i] for i in order]


def compute_genome_tree(genomes, hits, method="average", optimal_order=True,
                        zero_tolerance=0.0, verbose=True):
    """
    Returns (ordered_genomes, linkage_matrix, distance_matrix, peq_matrix,
             genome_names).

    IMPORTANT: distance_matrix, peq_matrix and genome_names are all in the
    ORIGINAL genome order. Only `ordered_genomes` is permuted. Pass
    `genome_names` (not the reordered names) to scipy's dendrogram as `labels`.
    """
    genome_names, peq_matrix = compute_peq_matrix(genomes, hits)
    distance_matrix = compute_distance_matrix(peq_matrix, zero_tolerance=zero_tolerance)

    if verbose:
        import pandas as pd

        df = pd.DataFrame(distance_matrix, index=genome_names, columns=genome_names)
        # %g rather than round(4): a distance of 2e-04 prints as 0.0000 under
        # rounding, which reads as "identical" while the tree still draws the
        # bracket for the residue that is actually there.
        print(df.to_string(float_format=lambda v: f"{v:.4g}"))

        off_diagonal = ~np.eye(len(genome_names), dtype=bool)
        print("PEQ  min/max (off-diagonal):",
              peq_matrix[off_diagonal].min(), peq_matrix[off_diagonal].max())
        print("Dist min/max (off-diagonal):",
              distance_matrix[off_diagonal].min(), distance_matrix[off_diagonal].max())

        # Report pairs that came out identical, near-identical or unrelated.
        # Near-identical pairs are the ones that draw as a small bracket rather
        # than a flat line; the printed value is what zero_tolerance must reach.
        for i in range(len(genome_names)):
            for j in range(i + 1, len(genome_names)):
                d_ij = distance_matrix[i, j]
                if d_ij == 0.0:
                    print(f"  identical: {genome_names[i]} vs {genome_names[j]}")
                elif d_ij < 0.01:
                    print(f"  near-identical: {genome_names[i]} vs "
                          f"{genome_names[j]} (distance {d_ij:.2e})")
                elif d_ij >= 1.0 - 1e-9:
                    print(f"  no shared genes: {genome_names[i]} vs {genome_names[j]}")

    linkage_matrix, order = build_tree(
        distance_matrix, method=method, optimal_order=optimal_order
    )
    ordered_genomes = order_genomes(genomes, order)

    return ordered_genomes, linkage_matrix, distance_matrix, peq_matrix, genome_names