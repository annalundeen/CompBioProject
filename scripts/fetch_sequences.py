import requests
from Bio import SeqIO
from io import StringIO
import sys
import os

input_file = sys.argv[1]
output_fasta = sys.argv[2]
output_motifs = sys.argv[3]

def parse_input(file):
    proteins = {}
    motifs = []
    section = None

    with open(file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("["):
                section = line
                continue

            if section == "[PROTEINS]":
                name, uid = [x.strip() for x in line.split(",")]
                proteins[name] = uid

            elif section == "[MOTIFS]":
                motifs.append(line)
    
    print("Proteins parsed:", proteins)
    print("Motifs parsed:", motifs)

    return proteins, motifs

def fetch_fasta(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    r = requests.get(url)

    if r.status_code != 200 or not r.text.startswith(">"):
        raise ValueError(f"Failed to fetch {uniprot_id}")

    return r.text

proteins, motifs = parse_input(input_file)

records = []

for name, uid in proteins.items():
    fasta = fetch_fasta(uid)
    record = SeqIO.read(StringIO(fasta), "fasta")
    record.id = name
    record.description = ""
    records.append(record)

os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
SeqIO.write(records, output_fasta, "fasta")

#save motifs separately too just in case
motif_path = os.path.join(os.path.dirname(output_fasta), "motifs.txt")

with open(output_motifs, "w") as f:
    for m in motifs:
        f.write(m + "\n")