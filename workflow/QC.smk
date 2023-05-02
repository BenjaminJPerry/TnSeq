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
        'output/00_QC/AllMultiQCReport.html',


rule fastqc_raw:
    input:
        fastq = 'fastq/{sample}.fastq.gz'
    output:
        html = 'output/00_QC/fastqc/{sample}_fastqc.html',
        zip = 'output/00_QC/fastqc/{sample}_fastqc.zip'
    conda:
        'fastqc'
        # 'docker://biocontainers/fastqc:v0.11.9_cv8'
    threads: 2
    message:
        'Running QC on reads: {wildcards.sample}\n'
    shell:
        'fastqc '
        '-o output/00_QC/fastqc/ '
        '-q '
        '-t {threads} '
        '{input.fastq}'


rule fastqc_bbduk:
    input:
        fastq = 'output/01_trimmed_reads/{sample}.bbduk.fastq.gz'
    output:
        html = 'output/00_QC/fastqc/{sample}.bbduk_fastqc.html',
        zip = 'output/00_QC/fastqc/{sample}.bbduk_fastqc.zip'
    conda:
        'fastqc'
        # 'docker://biocontainers/fastqc:v0.11.9_cv8'
    threads: 2
    message:
        'Running QC on reads: {wildcards.sample}\n'
    shell:
        'fastqc '
        '-o output/00_QC/fastqc_bbduk/ '
        '-q '
        '-t {threads} '
        '{input.fastq}'


rule fastqc_trimmed:
    input:
        fastq = 'output/01_trimmed_reads/{sample}.bbduk.trimmed.fastq.gz'
    output:
        html = 'output/00_QC/fastqc_trimmed/{sample}.bbduk.trimmed_fastqc.html',
        zip = 'output/00_QC/fastqc_trimmed/{sample}.bbduk.trimmed_fastqc.zip'
    conda:
        'fastqc'
        # 'docker://biocontainers/fastqc:v0.11.9_cv8'
    threads: 2
    message:
        'Running QC on reads: {wildcards.sample}\n'
    shell:
        'fastqc '
        '-o output/00_QC/fastqc_trimmed/ '
        '-q '
        '-t {threads} '
        '{input.fastq}'


rule multiQC:
    input:
        raw = expand('output/00_QC/fastqc/{sample}_fastqc.zip', sample = SAMPLES),
        bbduk = expand('output/00_QC/fastqc/{sample}.bbduk_fastqc.zip', sample = SAMPLES),
        trimming = expand('logs/cutadapt.{sample}.log', sample = SAMPLES),
        trimmed = expand('output/00_QC/fastqc_trimmed/{sample}.bbduk.trimmed_fastqc.zip', sample = SAMPLES),
        alignment = expand('logs/tntag_alignment.{sample}.log', sample = SAMPLES),
        filtering = expand('logs/pJG714_filter.{sample}.log', sample = SAMPLES),
        ecoli = expand('logs/ecoli_check.{sample}.log', sample = SAMPLES), 
    output:
        multiQC ='output/00_QC/AllMultiQCReport.html'
    conda:
        'multiqc'
        # 'docker://quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
    shell:
        'multiqc '
        '--profile-runtime '
        '-n output/00_QC/AllMultiQCReport '
        '-s '
        '-f '
        '--interactive '
        '{input.raw} {input.trimming} {input.trimmed} {input.alignment} {input.filtering} {input.ecoli} '

