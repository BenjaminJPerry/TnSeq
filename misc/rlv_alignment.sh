#!/bin/bash

echo "bowtie alignment and conversiont shell script..."
echo "This script will align a tntag fasta file to the Rlv3841 genome and then convert the alignment from .sam to .bed format."

cd "/usr/lib/tnseq/align/"

printf "Enter the tntag.fasta file for alignment: "
read TNTAGS

printf "Enter the output file for the .bed file: "
read OUTFILE


TEMP1="/tmp/alignment.sam"
touch $TEMP1
TEMP2="/tmp/alignment.bam"
touch $TEMP2

#align the reads using bowtie
echo "Aligning tntags to Rlv3841 reference genome."

echo `bowtie -f -n0 -l16 -m1 -p3 -S RLV "$TNTAGS" > "$TEMP1"`

echo "Alignment complete."

#convert the .sam alignment to .bam
echo "Converting .sam alignment into .bam format."

echo `samtools view -Sb "$TEMP1" > "$TEMP2"`

#convert the .bam to .bed format
echo "Converting .bam alignment into .bed format."

echo `bedtools bamtobed -i "$TEMP2" > "$OUTFILE"`

echo "Final alignment file printed to $OUTFILE."

cd "/home/benji"
