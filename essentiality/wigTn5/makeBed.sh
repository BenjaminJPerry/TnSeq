# 2018 Benjamin J Perry - Attribution-NonCommercial-ShareAlike 4.0 International
# (CC BY-NC-SA 4.0)
# Version: 1.2.0
# Maintainer: Benjamin J Perry
# Email: benjamin.perry@postgrad.otago.ac.nz
# Status: Functional

#!/bin/bash

# Overview:
# This shell script is designed to take a directory of single ended reads, each representing the data
# from a single Tn-Seq sequencing library and trim, align, and convert the data into useable format listed below.
#
# Output Files:
# 1. sorted .bam and .bai file for loading into IGV
# 2. bedfile of mapped reads for using in R scrip for strand separation and wig file generation
# 3. fastqc report of trimmed reads
# 4. log file of reads mapped to pJG714 and E. coli ST18 genome
#
# Synopsis or Pipeline Steps:
# 1. cutadapt is used to remove the Tn sequences and any polyC tails which may be on the 3' end.
# 2. bowtie2 aligns reads to pJG714 vector and save unaligned reads.
# 3. bowtie2 aligns reads to E. coli genome for observation of possible E. coli reads
# 4. Phix-pJG714 filter reads are aligned to the reference genome in question
# 5. Alignment output is converted, sorted, and indexed for visualization in IGV
# 6. Sorted .bam alignment file in converted into .bedfile for downstream analysis in R

#activate conda environment
clear
pwd
printf "\n\n\n"

#source activate py36

###Global variables needed for analysis
R7AWT=~/ref/R7A
R7ANS=~/Projects/TnSeq/ref/R7ANS
ECOLIREF=~/ref/ECOLI
PJG714REF=~/ref/PJG714
TA1CAT=/home/ronson/Projects/TnSeq/ref/TA1CAT
TA1CON=~/Projects/TnSeq/ref/TA1CON

###	Loop through each directory
for i in $(ls | cut -d '.' -f 1);
do

mkdir "$i"
mv "$i"*.gz "$i"

cd $i
# Working in Sample root dir
mkdir reads
mv "$i"*.gz reads
READS=reads/"$(ls reads)"

printf "\n\n\n"
printf "Trimming Reads...\n\n\n"
# Trim Tn5 ME and polyC; truncate reads at 50 bp
TRIMREADS=reads/"$i".trim.fastq.gz
cutadapt  -j 12 -g TGTGTATAAGAGACAG -l 50 -m 25 -e 0.2 --discard-untrimmed -o "$TRIMREADS" "$READS" > >(tee -i reads/"$i".cutadapt.log)

### Filter pJG714 reads from data
mkdir filter

JGSORTBAM=filter/pJG714."$i".sort.bam
JGSORTBAMBAI="$JGSORTBAM".bai
JGFILTREADS=reads/"$i".trim.filter.fastq.gz

printf "\n\n\n"
printf "Aligning Reads to pJG714 to Filter...\n\n\n"

bowtie2 -p 12 --fast -x "$PJG714REF" -U "$TRIMREADS" --un-gz "$JGFILTREADS" | samtools view -b | samtools sort -o "$JGSORTBAM"
samtools index "$JGSORTBAM" "$JGSORTBAMBAI"

### Align Reads to reference genome
mkdir alignment

SORTBAM=alignment/"$i".sort.bam
SORTBAMBAI="$SORTBAM".bai
UNALINREADS=reads/unaligned."$i".trim.filter.fastq.gz
BEDFILE=alignment/"$i".bed

printf "\n\n\n"
printf "Aligning Reads to R7A Genome...\n\n\n"

### Alignment to the desired reference genome
bowtie2 -p 12 --very-sensitive -x "$R7AWT" -U "$JGFILTREADS" --no-unal --un-gz "$UNALINREADS" | samtools view -q 30 -b | samtools sort -o "$SORTBAM"
samtools index "$SORTBAM" "$SORTBAMBAI"
bedtools bamtobed -i "$SORTBAM" > "$BEDFILE"

### Align remaining reads to E coli genome
mkdir contaminants

CONBAM=contaminants/ecoli."$i".sort.bam
CONBAMBAI="$CONBAM".bai
MYSTYREADS=contaminants/mystery.reads.fastq.gz

printf "\n\n\n"
printf "Aligning Remaining Reads to E. coli K12...\n\n\n"

bowtie2 -p 12 --fast -x "$ECOLIREF" -U "$UNALINREADS" --un-gz "$MYSTYREADS" | samtools view -b | samtools sort -o "$CONBAM"
samtools index "$CONBAM" "$CONBAMBAI"

#python /home/ronson/Projects/TnSeq/essentiality/tnScripts/wigScripts.py -F /home/ronson/Projects/TnSeq/ref/TA1.CAT.fasta -B "$BEDFILE" -O "$i".tn5.wig -Tn5

cd ..

done
