# 2019 Benjamin J Perry - Attribution-NonCommercial-ShareAlike 4.0 International
# (CC BY-NC-SA 4.0)
# Version: 0.1.0
# Maintainer: Benjamin J Perry
# Email: benjamin.perry@postgrad.otago.ac.nz
# Status: Dev
# Citation: TBD

# Description:
# 1. Process fastq files into .wig files for use with TRANSIT
# 2. Process fastq files into .wig files for viewing in IGV (Stranded and Total)
# 3. Process fastq files into .sort.sam file for analysis with TSAS


def


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser()

    #Check for blank call, return usage
    if len(sys.argv[1:]) == 0:
        parser.print_usage()
        parser.exit()

    ### Designate Inputs:

    # Input fastq file(s)
    parser.add_argument("-f", "--fastqFile", type=str, help="path to input fastq file")
    # Refence Genome Sequence
    parser.add_argument("-r", "--refenceFile", type=str, help="path to reference file")
    # Output File Name Base
    parser.add_argument("-o", "--outputDir", type=str, help="output directory")

    # Tn trim string
    parser.add_argument("-s", "--tnString", type=str, help="Tn substring to trim")

    # Threads
    parser.add_argument("-t", "--threads", type=str, help="processors to use")
    # Memory
    parser.add_argument("m", "--memory", type=str, help="Gb of memory to use")


    print('2018 Benjamin J Perry - (CC BY-NC-SA 4.0)')
    print('Version: 2.1.0')
    print('Email: benjamin.perry@postgrad.otago.ac.nz')
    print('Citation: TBD\n')
