# src/blast_links.py

# Definition of the BlastHit class. Using the __init__() built-in method to initiate the class. That way I can assign values to object properties.
    # Here I defined those parameters that can be found in a BLAST tabular output format 6:
        # qseqid   = query or source (gene) sequence id
        # sseqid   = subject or target (reference genome) sequence id
        # identity = percentage of identical positions
        # length   = alignment length (number of aligned amino acids)
        # qlen     = total length of the query sequence (in amino acids)
        # slen     = total length of the subject sequence (in amino acids)
        # qcoverage = fraction of the query sequence covered by the alignment (length / qlen)
        # scoverage = fraction of the subject sequence covered by the alignment (length / slen)
class BlastHit:
    def __init__(self, qseqid, sseqid, identity, length, qlen, slen):
        self.qseqid = qseqid
        self.sseqid = sseqid
        self.identity = identity
        self.length = length
        self.qlen = qlen
        self.slen = slen
        self.qcoverage = length / qlen  # % de la query alignée
        self.scoverage = length / slen  # % du subject aligné


# Definition of a fonction to extract blast data. 
def load_blast(blast_tsv, min_identity=0, min_coverage=0):
    """
    Load BLAST outfmt 6
    Return a list of all the blast hits
    """

    hits = []
    # Open all-vs-all blast tsv file
    with open(blast_tsv) as f:
        hits = []

        for i, line in enumerate(f):
            fields = line.strip().split("\t")
            
            # Debug first line
            if i == 0:
                print("Number of fields:", len(fields))
                print("First line:", fields)
            #Attribute each info to the appropriate variable
            qseqid = fields[0]
            sseqid = fields[1]
            pident = float(fields[2])
            length  = int(fields[3])
            qlen    = int(fields[4])
            slen    = int(fields[5])

            qcoverage = (length / qlen) * 100
            scoverage = (length / slen) * 100

            # For each line if the percentage of identity is lower than the decided threshold, 
            # this blast will be skip
            if pident < min_identity:
                continue
            if qcoverage < min_coverage:  # ← use the variable, not length/qlen
                continue
            
            #Clean up the gene name in the file
            qseqid = qseqid.removeprefix("Pharokka_")
            sseqid = sseqid.removeprefix("Pharokka_")

            if i == 0:
                print("After prefix removal:", qseqid, sseqid)
            
            #Remove the hit from a sequence with itself
            if qseqid == sseqid:
                continue
            
            #Add in the hits list the object BlastHit corresponding to a blast that add a higher identity 
            #score than the treshold and that it not a protein hit with itsself
            hits.append(BlastHit(qseqid, sseqid, pident, length, qlen, slen))


    #Return the list of object
    return hits
