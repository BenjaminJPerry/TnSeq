#!/bin/bash

cutadapt -q 20 --discard-untrimmed -e 0.1 -f fastq -a AGATCGGAAGAGCGTCGTGTAGGGAA *.fastq > *.3trim.fastq #Trimming the 3' Adpater

cutadapt --discard-untrimmed -m 14 -M 17 -e 0.1 -f fastq -g AGACCGGGGACTTATCATCCAACCTGT *.3trim.fastq > *.trim.fastq #Trimming the Mariner IR element

bowtie -n0 -l16 -m1 -p6 -S /home/benji/Desktop/Oxf_INSeq/Ref/RLV3841 *.trim.fastq > *.trim.sam #Aligning the trimmed reads to the RLV3841 Genome

samtools view -Sb *.trim.sam > *.trim.bam #Converting the .sam to a .bam

bedtools bamtobed -i *.trim.bam > *.trim.bed # Make a bedfile out of the .bam file

python /home/benji/Desktop/Oxf_INSeq/make_wig_exec.py *.trim.bed # Some python magic to convert the bedfile into 7 wig files

python /home/benji/Desktop/Oxf_INSeq/tn-hmm_1.03/tn-hmm.py -f AM236080.wig -gff /home/benji/Desktop/Oxf_INSeq/Ref/AM236080.gff3 > AM236080.hmm.out # Running the HMM

python /home/benji/Desktop/Oxf_INSeq/tn-hmm_1.03/process_genes.py -f AM236080.hmm.out > AM236080.process.txt # Processing the HMM with the updated NCBI annotation


python /home/benji/Desktop/Oxf_INSeq/tn-hmm_1.03/tn-hmm.py -f AM236086.wig -gff /home/benji/Desktop/Oxf_INSeq/Ref/AM236086.gff3 > AM236086.hmm.out

python /home/benji/Desktop/Oxf_INSeq/tn-hmm_1.03/process_genes.py -f AM236086.hmm.out > AM236086.process.txt


python /home/benji/Desktop/Oxf_INSeq/tn-hmm_1.03/tn-hmm.py -f AM236085.wig -gff /home/benji/Desktop/Oxf_INSeq/Ref/AM236085.gff3 > AM236085.hmm.out

python /home/benji/Desktop/Oxf_INSeq/tn-hmm_1.03/process_genes.py -f AM236085.hmm.out > AM236085.process.txt


python /home/benji/Desktop/Oxf_INSeq/tn-hmm_1.03/tn-hmm.py -f AM236084.wig -gff /home/benji/Desktop/Oxf_INSeq/Ref/AM236084.gff3 > AM236084.hmm.out

python /home/benji/Desktop/Oxf_INSeq/tn-hmm_1.03/process_genes.py -f AM236084.hmm.out > AM236084.process.txt


python /home/benji/Desktop/Oxf_INSeq/tn-hmm_1.03/tn-hmm.py -f AM236083.wig -gff /home/benji/Desktop/Oxf_INSeq/Ref/AM236083.gff3 > AM236083.hmm.out

python /home/benji/Desktop/Oxf_INSeq/tn-hmm_1.03/process_genes.py -f AM236083.hmm.out > AM236083.process.txt


python /home/benji/Desktop/Oxf_INSeq/tn-hmm_1.03/tn-hmm.py -f AM236082.wig -gff /home/benji/Desktop/Oxf_INSeq/Ref/AM236082.gff3 > AM236082.hmm.out

python /home/benji/Desktop/Oxf_INSeq/tn-hmm_1.03/process_genes.py -f AM236082.hmm.out > AM236082.process.txt


python /home/benji/Desktop/Oxf_INSeq/tn-hmm_1.03/tn-hmm.py -f AM236081.wig -gff /home/benji/Desktop/Oxf_INSeq/Ref/AM236081.gff3 > AM236081.hmm.out

python /home/benji/Desktop/Oxf_INSeq/tn-hmm_1.03/process_genes.py -f AM236081.hmm.out > AM236081.process.txt

sleep 10 # Let the shell catch up to the output files before moving on to the next directory
