# 2023 Benjamin J Perry
# MIT License
# Copyright (c) 2023 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

# configfile: "config/config.yaml"


import os

wildcard_constraints:
   sample = '\w+',


(SAMPLES,) = glob_wildcards("output/04_aligned_beds/{sample}.bed")

(REPLICONS,) = glob_wildcards("ref/embl/{replicon}.embl")


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
        expand("output/07_tradis_output/{sample}.{replicon}.tradis_gene_insert_sites.tsv.essen.csv", sample = SAMPLES, replicon = REPLICONS),

rule gene_insertions:
    input:
        annotation="ref/embl/{replicon}.embl",
        insertions="output/05_tradis_plots/{sample}.{replicon}.insert_site_plot.gz"
    output:
        tradis_gene_counts="output/07_tradis_output/{sample}.{replicon}.tradis_gene_insert_sites.tsv"
    log:
        "logs/tradis_gene_insert_sites.{sample}.{replicon}.log"
    threads: 2
    conda:
        "tradis"
    shell:
        "tradis_gene_insert_sites "
        #"-o {wildcards.sample}.{wildcards.replicon}.tradis_gene_insert_sites.csv "
        "-trim5 0.1 "
        "-trim3 0.1 "
        "{input.annotation} "
        "{input.insertions} "
        "&& " # JANKY because tradis made me do it...
        "mv {wildcards.sample}.{wildcards.replicon}.tradis_gene_insert_sites.csv {output.tradis_gene_counts} "

rule gene_essentiality:
    input:
        tradis_gene_counts="output/07_tradis_output/{sample}.{replicon}.tradis_gene_insert_sites.tsv",
    output:
        states = "output/07_tradis_output/{sample}.{replicon}.tradis_gene_insert_sites.tsv.all.csv",
        ambiguous = "output/07_tradis_output/{sample}.{replicon}.tradis_gene_insert_sites.tsv.ambig.csv",
        essential = "output/07_tradis_output/{sample}.{replicon}.tradis_gene_insert_sites.tsv.essen.csv",
        change_point_plot = "output/07_tradis_output/{sample}.{replicon}.tradis_gene_insert_sites.tsv.QC_and_changepoint_plots.pdf",
    log:
        "logs/tradis_gene_insert_sites.{sample}.{replicon}.log"
    threads: 2
    conda:
        "tradis"
    shell:
        "tradis_essentiality.R {input.tradis_gene_counts}"
