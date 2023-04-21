# 2023 Benjamin J Perry
# MIT License
# Copyright (c) 2023 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

# configfile: "config/config.yaml"


import os


(FIDs,) = glob_wildcards("results/02_kneaddata/{sample}.fastq")


onstart:
    print(f"Working directory: {os.getcwd()}")
    print("TOOLS: ")
    os.system('echo "  bash: $(which bash)"')
    os.system('echo "  PYTHON: $(which python)"')
    os.system('echo "  CONDA: $(which conda)"')
    os.system('echo "  SNAKEMAKE: $(which snakemake)"')
    print(f"Env TMPDIR = {os.environ.get('TMPDIR', '<n/a>')}")
    os.system('echo "  PYTHON VERSION: $(python --version)"')
    os.system('echo "  CONDA VERSION: $(conda --version)"')


rule all:
    input:
        "results/centrifuge.counts.tsv",
        "results/centrifuge.counts.biom",

        "results/kraken2.counts.tsv",
        "results/kraken2.counts.biom",

        "results/bracken.k2.counts.tsv",
        "results/bracken.k2.counts.biom",

        # expand("results/03_humann3Uniref50EC/{sample}_pathcoverage.tsv", sample=FIDs),


localrules:
    generateCentrifugeSampleSheet,


rule generateCentrifugeSampleSheet:
    output:
        sampleSheet = "resources/centrifugeSampleSheet.tsv",
    threads: 2
    shell:
        "./workflow/scripts/generate_centrifuge_sample_sheet.sh -d results/02_kneaddata -p fastq -o {output.sampleSheet} "


rule centrifugeGTDB:
    input:
        sampleSheet = "resources/centrifugeSampleSheet.tsv",
    output:
        out = expand("results/03_centrifuge/{sample}.GTDB.centrifuge", sample = FIDs),
        report = expand("results/03_centrifuge/{sample}.GTDB.centrifuge.report", sample = FIDs),
    log:
        "logs/centrifuge.GTDB.multi.log",
    conda:
        "centrifuge"
    threads: 32
    resources:
        mem_gb = lambda wildacards, attempt: 140 + ((attempt - 1) + 20),
        time = "06:00:00",
    shell:
        "centrifuge "
        "-x /bifo/scratch/2022-BJP-GTDB/2022-BJP-GTDB/centrifuge/GTDB "
        "--sample-sheet {input.sampleSheet} "
        "-t "
        "--threads {threads} "
        "&> {log} "
