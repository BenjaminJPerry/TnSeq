# 2023 Benjamin J Perry
# MIT License
# Copyright (c) 2023 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: "config/config.yaml"


import os


(SAMPLES,) = glob_wildcards("fastq/{sample}.fastq.gz")


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
        expand("04_aligned_beds/{sample}.bed", sample=SAMPLES),


rule cutadapt:
    input:
        fastq="fastq/{sample}.fastq.gz"
    output:
        trimmed="01_trimmed_reads/{sample}.trimmed.fastq.gz"
    log:
        "logs/cutadapt.{sample}.log"
    threads: 12
    conda:
        "workflow/envs/cutadapt.yaml"
    params:
        IR=config["IR_ELEMENT"]
    shell:
        "cutadapt "
        "-j {threads} "
        "-g {params.IR} "
        "-l 50 "
        "-m 25 "
        "-e 0.2 "
        "--discard-untrimmed "
        "-o {output.trimmed} "
        "{input.fastq} "
        "2>&1 | tee {log}"


rule pJG714_Filter:
    input:
        trimmed="01_trimmed_reads/{sample}.trimmed.fastq.gz"
    output:
        pJG714="02_pJG714_filtering/{sample}.pJG714.fastq.gz",
        filtered="02_pJG714_filtering/{sample}.trimmed.filtered.fastq.gz"
    log:
        "logs/pJG714_filter.{sample}.log"
    threads: 8
    conda:
        "workflow/envs/bowtie2.yaml"
    params:
        PJG714=config["BOWTIE_PJG714REF"]
    shell:
        "bowtie2 "
        "-p {threads} "
        " --fast "
        "-x {params.PJG714} "
        "-U {input.trimmed} "
        "--un-gz {output.filtered} "
        "--al-gz {output.pJG714} "
        "> /dev/null"


rule align_tags:
    input:
        filtered="02_pJG714_filtering/{sample}.trimmed.filtered.fastq.gz",
    output:
        aligned_bam="03_aligned_bams/{sample}.trimmed.filtered.sorted.bam",
        aligned_bai="03_aligned_bams/{sample}.trimmed.filtered.sorted.bam.bai",
        unaligned="03_aligned_bams/{sample}.unaligned.fastq.gz",
    log:
        "logs/tntag_alignment.{sample}.log"
    threads: 8
    conda:
        "workflow/envs/bowtie2.yaml"
    params:
        REFERENCE=config["BOWTIE_REFERENCE"]
    shell:
        "bowtie2 "
        "-p {threads} "
        "--sensitive "
        "-x {params.REFERENCE} "
        "-U {input.filtered} "
        "--no-unal "
        "--un-gz {output.unaligned} | "
        "samtools "
        "sort "
        "--write-index"
        "-o {output.aligned_bam} "
        "2>&1 | tee {log}"

rule bam_to_bed:
    input:
        aligned_bam="03_aligned_bams/{sample}.trimmed.filtered.sorted.bam",
    output:
        aligned_bed="04_aligned_beds/{sample}.bed"
    log:
        "logs/bam_to_bed.{sample}.log"
    threads: 8
    conda:
        "workflow/envs/bedtools.yaml"
    shell:
        "bedtools "
        "bamtobed "
        "-i {input.aligned_bam} > {output.aligned_bed}"

