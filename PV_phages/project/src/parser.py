# src/parser.py
from Bio import SeqIO

# Define Genome class. Definition of the BlastHit class. Using the __init()__ built-in method to initiate the class. That way I can assign values to object properties. 
# That class will contain the information about the genome.
    #Name: Genome's name
    #Lenght: Genome's lenght
    #Genes: Number of genes
class Genome:
    def __init__(self, name, length, genes):
        self.name = name
        self.length = length
        self.genes = genes

# Define the Gene class. Using the __init()__ built-in method to initiate the class. That way I can assign values to object properties. 
# That class will contain the information about the gene.
    #Start: first coordonate of the gene
    #End : Last coordonate of the gene 
    #locus: Locus on which the gene is 
class Gene:
    def __init__(self, start, end, strand, locus, function):
        self.start = start
        self.end = end
        self.strand = strand
        self.locus_tag = locus
        self.function = function    
        # Return a printable representation of the object
    def __repr__(self):
        return f"Gene({self.locus_tag}, {self.start}-{self.end}, strand={self.strand})"
       
#Define a function 
def parse_genbank(gbk_file):
    """
    
    """
    #Create an empty list. That list will contain all the genomes 
    genomes = []
    # Iterate over all the gbk file 
    for record in SeqIO.parse(gbk_file, "genbank"):
        # genes = []
        #Each gbk file will be attribute to a genome object. Using the Seqrecord object created with SeqIO.parse, we are able the recovert 
        #the name, the lenght and a list that will contain all the genes coordonates
        genome = Genome(
            name=record.name,
            length=len(record.seq),
            genes=[]
        )
        # Iterate over all the features which is another list property, and it contains SeqFeature objects
        for feature in record.features:
        
            #For all CDS the gene coordinate will be extract and put in the gene list, as well as the strand on which the gene is
            if feature.type == "CDS":
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand
                function = feature.qualifiers.get("function", ["unknown"])[0]
                #Search in he qualifiers directory that contains all the annotations for each CDS for the annotation 'locus_tag'. If it doesn't exist
                #the code return 'NA. Since qualifier stores list, I force
                # the search and the replacement to the first element.
                # locus = feature.qualifiers.get("locus_tag", ["NA"])[0]


                locus = feature.qualifiers.get("locus_tag", [None])[0]

                if locus is None:
                    locus = f"{record.name}_{start}_{end}"
                locus = locus.removeprefix("Pharokka_")
                # Defying the gene object for each gene (because it is in the loop)
                gene = Gene(start, end, strand, locus, function)

                #Adding the genome property to the gene object and assigning the genome object itself into this new property
                gene.genome = genome        
                
                #Adding every gene to the genes list already created
                genome.genes.append(gene)

        # Add every genome to the genomes list already created   
        genomes.append(genome)

    return genomes
