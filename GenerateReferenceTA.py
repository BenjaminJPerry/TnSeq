#!/usr/bin/env python

import os
import sys
import argparse

def wigRefGen(fastafile, outfile):
    '''
    Reads in a single fasta sequence from a file. Identifies the position of the T in ever TA motif.
    Returns a wig files with every TA site and the count column set to 100.
    Fasta header line is used to the wig config line.

    :param fastafile: input fasta file.
    :param outfile: outfile to write reference 'TA' wig track to.
    :return: null
    '''


if __name__ == '__main__':
    # Parse arguments

    # Read input fasta file; check for multiple header lines.

    # Open the output file: outFile

    # Compute Wig Reference Track
    wigRefGen(inFile, outFile)
