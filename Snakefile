#final target rule
rule all:
    input: #final output files
        "results/blast/msa.tsv", #output from blast
        "results/meme/meme_clean.tsv", #ouput from meme
        "results/final/simple_unique_hits.tsv", #minimal motif hits
        "results/final/ranked_results.tsv" #motif hits with metrics

#fetch sequences and motifs 
rule fetch:
    input: #config file with proteins and motifs
        "config/input.txt"
    output:
        fasta="data/raw/sequences.fasta", #output fastas for proteins
        motifs="data/raw/motifs.txt" #output motifs of interest
    shell: #run python script to parse inout and get sequences
        "python scripts/fetch_sequences.py {input} {output.fasta} {output.motifs}"

#run blast motifs
rule blast:
    input: #uses sequences and motifs
        fasta="data/raw/sequences.fasta",
        motifs="data/raw/motifs.txt"
    output: #output aligned blast results
        "results/blast/msa.tsv"
    shell: #runs blast pipeline and produces tsv of results
        "python scripts/run_blast.py {input.fasta} {input.motifs} results/blast"

#run meme motifs
rule meme_raw:
    input:
        fasta="data/raw/sequences.fasta", #sequences
        motifs="data/raw/motifs.txt", #motifs
        msa="results/blast/msa.tsv" #blast results
    output: #output filtered meme results
        "results/meme/meme_clean.tsv"
    shell: #run meme workflow and output results table
        "python scripts/run_full_meme.py {input.fasta} {input.motifs} {input.msa} {output}"

#combine blast and meme results
rule combine:
    input: #blast and meme output
        blast="results/blast/msa.tsv",
        meme="results/meme/meme_clean.tsv"
    output:
        "results/final/simple_unique_hits.tsv", #simple results
        "results/final/ranked_results.tsv" #expanded results
    shell: #script to merge results from methods, filter results, output scored results
        "python scripts/combine_results.py {input.blast} {input.meme} results/final"

#cleanup rule
rule clean:
    shell: #remove all generated result directories and raw intermediate files
        """
        rm -rf results/blast results/meme results/final
        rm -f data/raw/sequences.fasta data/raw/motifs.txt
        """