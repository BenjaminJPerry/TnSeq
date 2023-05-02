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
        expand("output/05_tradis_plots/.{sample}", sample=SAMPLES),
        expand("output/00_ECOLI_check/{sample}.Ecoli.K12.bam", sample=SAMPLES)


rule cutadapt:
    input:
        reads = 'fastq/{sample}.fastq.gz',
    output:
        trimmed="output/01_trimmed_reads/{sample}.trimmed.fastq.gz"
    log:
        "logs/cutadapt.{sample}.log"
    threads: 12
    conda:
        "envs/cutadapt.yaml"
    params:
        IR=config["IR_ELEMENT"]
    shell:
        "cutadapt "
        "-j {threads} "
        "-g {params.IR} " #5' clip
        "-O 12 " # minimum overlap
        #"-u 9 "
        #"-l 50 "
        #"-m 25 "
        #"-q 15 "
        "-e 0.2 " # missmatching allowed
        #"--discard-untrimmed "
        "-o {output.trimmed} "
        "{input.reads} "
        "2>&1 | tee {log}"



rule bbduk:
    input:
        reads="output/01_trimmed_reads/{sample}.trimmed.fastq.gz" 
    output:
        bbdukReads = 'output/01_trimmed_reads/{sample}.trimmed.bbduk.fastq.gz'
    log:
        'logs/bbduk.{sample}.log'
    conda:
        'bbduk'
    threads:4
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in={input.reads} '
        # 'entropy=0.3 '
        # 'entropywindow=50 '
        'trimpolygright=5 '
        'qtrim=r '
        'trimq=15 '
        'out={output.bbdukReads} '
        '2>&1 | tee {log}'


rule pJG714_Filter:
    input:
        trimmed="output/01_trimmed_reads/{sample}.trimmed.bbduk.fastq.gz"
    output:
        pJG714="output/02_pJG714_filtering/{sample}.pJG714.fastq.gz",
        filtered="output/02_pJG714_filtering/{sample}.trimmed.bbduk.filtered.fastq.gz"
    log:
        "logs/pJG714_filter.{sample}.log"
    threads: 8
    conda:
        "envs/bowtie2.yaml"
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
        "> /dev/null "
        "2> {log}"


rule align_tags:
    input:
        filtered="output/02_pJG714_filtering/{sample}.trimmed.bbduk.filtered.fastq.gz",
    output:
        aligned_bam="output/03_aligned_bams/{sample}.trimmed.bbduk.filtered.sorted.bam",
        aligned_bai="output/03_aligned_bams/{sample}.trimmed.bbduk.filtered.sorted.bam.bai",
        unaligned="output/03_aligned_bams/{sample}.unaligned.fastq.gz",
    log:
        "logs/tntag_alignment.{sample}.log"
    threads: 8
    conda:
        "envs/bowtie2.yaml"
    params:
        REFERENCE=config["BOWTIE_REFERENCE"]
    shell:
        "bowtie2 "
        "-p {threads} "
        "--sensitive "
        "-x {params.REFERENCE} "
        "-U {input.filtered} "
        "--no-unal "
        "--un-gz {output.unaligned} "
        "2> {log} "
        "| "
        "samtools "
        "sort "
        "-o {output.aligned_bam} "
        "&& "
        "samtools index {output.aligned_bam} "


rule ECOLI_Check:
    input:
        leftovers="output/03_aligned_bams/{sample}.unaligned.fastq.gz"
    output:
        temp("output/00_ECOLI_check/{sample}.Ecoli.K12.bam"),
        unaligned = "output/00_ECOLI_check/{sample}.Ecoli.K12.unaligned.fastq.gz"
    log:
        "logs/ecoli_check.{sample}.log"
    threads: 8
    conda:
        "envs/bowtie2.yaml"
    params:
        ECOLI=config["BOWTIE_ECOLIREF"]
    shell:
        "bowtie2 "
        "-p {threads} "
        " --fast "
        "-x {params.ECOLI} "
        "-U {input.leftovers} "
        "--un-gz {output.unaligned} "
        "> {output} "
        "2> {log}"


rule bam_to_bed:
    input:
        aligned_bam="output/03_aligned_bams/{sample}.trimmed.bbduk.filtered.sorted.bam",
    output:
        aligned_bed="output/04_aligned_beds/{sample}.bed"
    log:
        "logs/bam_to_bed.{sample}.log"
    threads: 8
    conda:
        "envs/bedtools.yaml"
    shell:
        "bedtools "
        "bamtobed "
        "-i {input.aligned_bam} > {output.aligned_bed}"


rule bed_to_insertionplot:
    input:
        aligned_bam="output/03_aligned_bams/{sample}.trimmed.bbduk.filtered.sorted.bam",
        aligned_bed="output/04_aligned_beds/{sample}.bed"
    output:
        tradis_insertionplot_semaphore="output/05_tradis_plots/.{sample}"
    log:
        "logs/bed_to_insertionplot.{sample}.log"
    threads: 2
    conda:
        "envs/insertionPlots.yaml"
    shell:
        "Rscript workflow/scripts/makeInsertionplots.R {input.aligned_bam} {input.aligned_bed} {wildcards.sample} "
        "> {log} "
        "&& touch {output.tradis_insertionplot_semaphore} "


