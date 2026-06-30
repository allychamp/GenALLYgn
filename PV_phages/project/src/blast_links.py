# src/blast_links.py

# Definition of the BlastHit class. Using the __init()__ built-in method to initiate the class. That way I can assign values to object properties. 
    #Here I defiened those parameters that can be found in a Blastn tabular output format 6: 
        #qseqid = query or source (gene) sequence id
        #sseqid = subject or target (reference genome) sequence id
        #identity =  percentage of identical positions
class BlastHit:
    def __init__(self, qseqid, sseqid, identity):
        self.qseqid = qseqid
        self.sseqid = sseqid
        self.identity = identity

# Definition of a fonction to extract blast data. 
def load_blast(blast_tsv, min_identity=0):
    """
    Load BLAST outfmt 6
    Return a list of all the blast hits
    """

    hits = []
    # Open all-vs-all blast tsv file
    with open(blast_tsv) as f:
        #Iterate over each line 
        for line in f:
            # if line.startswith("#") or not line.strip():
            #     continue

            #Separate each field using the tab separator in the file
            fields = line.strip().split("\t")

            #Attribute each info to the appropriate variable
            qseqid = fields[0]
            sseqid = fields[1]
            pident = float(fields[2])

            # For each line if the percentage of identity is lower than the decided threshold, 
            # this blast will be skip
            if pident < min_identity:
                continue

            #Clean up the gene name in the file
            qseqid = qseqid.removeprefix("Pharokka_")
            sseqid = sseqid.removeprefix("Pharokka_")

            #Remove the hit from a sequence with itself
            if qseqid == sseqid:
                continue
            
            #Add in the hits list the object BlastHit corresponding to a blast that add a higher identity 
            #score than the treshold and that it not a protein hit with itsself
            hits.append(BlastHit(qseqid, sseqid, pident))

    #Return the list of object
    return hits
