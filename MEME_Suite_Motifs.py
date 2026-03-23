###############################
# Used the Module Pymemesuite which allows you to access meme suit's tool FIMO:Find Individual Motif Occurrences
    #FIMO: scans a set of sequences for individual matches to motifs you provide 
################################
#LINK TO SITE WITH INSTALLION INSTRUCTIONS FOR MODULE
#https://pypi.org/project/pymemesuite/#:~:text=%F0%9F%94%A7%20Installing,qvalue%20)
# from command line installed the module pymemesuite with the command:
#pip install pymemesuite 

import requests 
import Bio.SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO

####################
# GET MUSCLE PROTEIN SEQUENCES
#####################
proteins = {
    "ABLIM1": "O14639",
    "MYBPC": "Q14896",
    "MYL2": "P10916"
}

def fetch_fasta(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    return requests.get(url).text

sequences = {}
fasta_list = []

for name, uid in proteins.items():
    fasta = fetch_fasta(uid)
    fasta_list.append(fasta)
    record = SeqIO.read(StringIO(fasta), "fasta")
    sequences[name] = str(record.seq)
    #seq_list.append(sequences[name])

fasta_files = "/n".join(fasta_list)
fasta_input = fasta_files.replace("/n", "")

with open ("muscle_seq.fna", "w") as f: 
    f.write(fasta_input)

#####################
#GET INPUT BACKGROUND MATRIX 
#####################
from Bio import SeqIO
'''from Bio.Align import substitution_matrices'''
#import numpy as np

# function format input for motifs
# edit to include background frequency lines
def write_motif_file(motifs, out_file):
    alphabet = "ACDEFGHIKLMNPQRSTVWY"

    with open(out_file, "w") as f:
    # write header
        f.write("MEME version 5\n")
        f.write(f"ALPHABET= {alphabet}\n\n")

        for i, seq in enumerate(motifs):    
        # motif identifier and alternate name
            f.write(f"MOTIF {seq} motif_{i+1}\n")

        # letter-probability matrix header
            f.write(f"letter-probability matrix: alength= {len(alphabet)} w= {len(seq)} nsites= 1\n")

        # design probability matrix
        # shows the probabilities of each amino acid at each position in motif
        # how to include scoring / chemical group similarity? this would modify code to account for weak motifs
            for i in seq:
                row = []
                for aa in alphabet:
                    if i == aa:
                        row.append("1.0000")
                    else:   
                        row.append("0.0000")
                f.write(" " + " ".join(row) + "\n")

motifs = ["GLALSDLIQKYFF"]
#"LSDLIQ", "LSSLIQ" other motifs to check once known to work fully 
write_motif_file(motifs, "motifs_input.txt")

############################
#START OF MEME SUITE MODULE:FIMO MOTIF SEARCH 
############################


#GENERATING MOTIF FREQUENCY MATRIX# 
from pymemesuite.common import MotifFile
with MotifFile("motifs_input.txt") as motif_file:
    motif = motif_file.read()

print(motif.name.decode())
print(motif.consensus)

for row in motif.frequencies:
    print(" ".join(f'{freq:.2f}' for freq in row))

#GENERATING MOTIF MATCHES#
from pymemesuite.common import Sequence 
from pymemesuite.fimo import FIMO 

sequences = [
    Sequence(str(record.seq), name=record.id.encode())
    for record in Bio.SeqIO.parse("muscle_seq.fna", "fasta")
]

fimo = FIMO(both_strands=False)
pattern = fimo.score_motif(motif,sequences, motif_file.background)

for m in pattern.matched_elements: 
    print(
    m.source.accession.decode(),
    m.start,
    m.stop,
    m.strand,
    m.score,
    m.pvalue,
    m.qvalue
    )

pos_1 = m.start 
pos_2 = m.stop 
match_name = m.source.accession.decode()

seq_id_dict = {}
for record in Bio.SeqIO.parse("muscle_seq.fna", "fasta"):
    seqs = str(record.seq)
    name = record.id
    seq_id_dict.update({name:seqs})

match_seq = seq_id_dict[match_name]
print(match_seq[pos_1:pos_2])