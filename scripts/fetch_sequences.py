#IMPORT PACKAGES
import requests #get sequences from uniprot
from Bio import SeqIO #handle sequence files
from io import StringIO #treat strings like file objects
import sys #command line
import os #paths/directories

#INPUTS
input_file = sys.argv[1]
output_fasta = sys.argv[2]
output_motifs = sys.argv[3]

#FUNCTION: PARSE INPUT
def parse_input(file):
    proteins = {} #dictionary with protein name and uniprot id
    motifs = [] #list of motif sequences
    section = None #track section of file

    with open(file) as f: #open config file and read
        for line in f:
            line = line.strip() #remove whitespace
            if not line:
                continue #skip empty lines
            if line.startswith("["): #find section headers
                section = line
                continue
            #parse protein lines
            if section == "[PROTEINS]":
                name, uid = [x.strip() for x in line.split(",")]
                proteins[name] = uid
            #parse motif lines
            elif section == "[MOTIFS]":
                motifs.append(line)
    #check if it worked
    print("Proteins parsed:", proteins)
    print("Motifs parsed:", motifs)

    return proteins, motifs

#FUNCTION: fetch fasta
def fetch_fasta(uniprot_id):
    #create url for getting fasta from uniprot
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    r = requests.get(url) #send get request
    #check for failure
    if r.status_code != 200 or not r.text.startswith(">"):
        raise ValueError(f"Failed to fetch {uniprot_id}")

    return r.text #get fasta text

#FETCH DATA
proteins, motifs = parse_input(input_file) #parse input file into proteins and motifs

records = [] #store seqrecord objects

#loop over each protein and get sequence
for name, uid in proteins.items():
    fasta = fetch_fasta(uid) #get fasta string from uniprot
    record = SeqIO.read(StringIO(fasta), "fasta") #convert string into seqrecord object
    record.id = name #replace id with readable protein name
    record.description = "" #remove long description
    records.append(record) #store record

#make sure output directory exists
os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
#write all sequences to fasta file
SeqIO.write(records, output_fasta, "fasta")

#save motifs separately too just in case
motif_path = os.path.join(os.path.dirname(output_fasta), "motifs.txt")

#write output
with open(output_motifs, "w") as f:
    for m in motifs:
        f.write(m + "\n")