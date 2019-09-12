#!/bin/bash

printf "Enter the path to a .fastq file -->"
read FASTQIN
echo "$FASTQIN read..."

TEMP2="/tmp/adapt_trim.fastq"
touch $TEMP2

echo `cutadapt -q 20 --discard-untrimmed -e 0.1 -f fastq -a AGATCGGAAGAGCGTCGTGTAGGGAA $FASTQIN > $TEMP2`

#Trim the IR element Sequence from the 5' end of each read

printf "Enter the path to a desired ouput file -->"
read TRIM_OUT
echo "$TRIM_OUT selected as output file."

echo `cutadapt --discard-untrimmed -m 14 -M 17 -e 0.1 -f fastq -g AGACCGGGGACTTATCATCCAACCTGT "$TEMP2" > "$TRIM_OUT"`
