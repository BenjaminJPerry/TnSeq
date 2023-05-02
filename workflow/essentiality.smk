# 2023 Benjamin J Perry
# MIT License
# Copyright (c) 2023 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

# configfile: "config/config.yaml"


import os


(SAMPLES,) = glob_wildcards("output/04_aligned_beds/{sample}.bed")
(REPLICONS,) = glob_wildcards("ref/{replicon}.embl")


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
        "",

rule gene_insertion_table:
    input:
        aligned_bam="output/03_aligned_bams/{sample}.trimmed.filtered.sorted.bam",
        aligned_bed="output/04_aligned_beds/{sample}.bed"
    output:
        tradis_insertionplot_semaphore="output/05_tradis_plots/.{sample}"
    log:
        "logs/bed_to_insertionplot.{sample}.log"
    threads: 2
    conda:
        "envs/insertionPlots.yaml"
    shell:
        "Rscript workflow/scripts/makeInsertionplots.R {input.aligned_bam} {input.aligned_bed} {wildcards.sample} 2>&1 {log} "
        "&& touch {output.tradis_insertionplot_semaphore} "
