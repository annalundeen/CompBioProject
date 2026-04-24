rule all:
    input:
        "results/blast/msa.tsv",
        "results/meme/meme_clean.tsv",
        "results/final/simple_unique_hits.tsv",
        "results/final/ranked_results.tsv"

rule fetch:
    input:
        "config/input.txt"
    output:
        fasta="data/raw/sequences.fasta",
        motifs="data/raw/motifs.txt"
    shell:
        "python scripts/fetch_sequences.py {input} {output.fasta} {output.motifs}"

rule blast:
    input:
        fasta="data/raw/sequences.fasta",
        motifs="data/raw/motifs.txt"
    output:
        "results/blast/msa.tsv"
    shell:
        "python scripts/run_blast.py {input.fasta} {input.motifs} results/blast"

rule meme_raw:
    input:
        fasta="data/raw/sequences.fasta",
        motifs="data/raw/motifs.txt",
        msa="results/blast/msa.tsv"
    output:
        "results/meme/meme_clean.tsv"
    shell:
        "python scripts/run_full_meme.py {input.fasta} {input.motifs} {input.msa} {output}"

rule combine:
    input:
        blast="results/blast/msa.tsv",
        meme="results/meme/meme_clean.tsv"
    output:
        "results/final/simple_unique_hits.tsv",
        "results/final/ranked_results.tsv"
    shell:
        "python scripts/combine_results.py {input.blast} {input.meme} results/final"

rule clean:
    shell:
        """
        rm -rf results/blast results/meme results/final
        rm -f data/raw/sequences.fasta data/raw/motifs.txt
        """