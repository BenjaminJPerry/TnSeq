# 2020 Benjamin J Perry - Attribution-NonCommercial-ShareAlike 4.0 International
# (CC BY-NC-SA 4.0)
# Version: 2.3.0
# Maintainer: Benjamin J Perry
# Email: benjamin.perry@postgrad.otago.ac.nz
# Status: Dev
# Citation: TBD

# Description:
# Python module for handling Tn5 related Tn-seq data.

def tn5RefWig(fastafile):
    """(reference.fasta) -> tn5Wig[pos][count = 0]
    :param fastafile: String; Input reference fasta file.
    :return: List[][]; [position][count], with count = 0.

    Reads in a single fasta sequence from a file.
    Returns a wig file of every position with [count] = 0.
    Fasta header line is used in wig config line.
    """

    print('Computing Reference Tn5 .wig file for: ' + fastafile + '\n')
    inputFile = open(fastafile, 'r')
    fastaData = inputFile.read()

    # Check fastafile for a single header line.
    if fastaData.count('>') > 1:
        print('Error: More than one fasta header line in ' + fastafile + '\n')
        raise SystemExit(1)

    # Process the fasta file into a header and a sequence variables
    fastaData = fastaData.splitlines()
    headerLine = fastaData.pop(0)
    headerLine = headerLine.split()[0]
    headerLine = headerLine[1:]
    fastaSequence = ''.join(fastaData)

    printHeader = ['variableStep chrom=' + str(headerLine)]
    Wig = []
    Wig.append(printHeader)
    i = 1

    # Make a list containing each position in genome
    for char in fastaSequence:
        entry = [i, 0]
        Wig.append(entry)
        i += 1

    return Wig


def tn5BedfileToTntags(bedfile):
    """ ('path/to/bedfile.bed') -> tntagsList[]
    :param bedfile: String; path to an input bedfile.
    :return: List[pos]; Tn5 insertion positions of all reads.

    Parses a bedfile from Tn5 Tnseq and returns a list of
    genomic coordinates offset by 9 bp for  Tn5 insertion positions.
    Takes into account the strand the read was mapped to.
    """
    # Open input bedfile
    bedfileIn = open(bedfile, 'r')

    # Process input bedfile
    bedfileData = []
    for line in bedfileIn:
        items = line.rstrip('\r\n').split('\t')
        items = [item.strip() for item in items]
        bedfileData.append(items)

    # Process tntag start positions
    tntagsList = []
    i = 0
    errors = 0
    for line in bedfileData:
        if str(line[5]) == '+':
            position = int(line[1])
            position = position + 1  # offsets python counting from 0
            position = position + 9  # offsets 9 bp duplication from Tn5
            tntagsList.append(position)
            i += 1
        elif str(line[5]) == '-':
            position = int(line[2])
            position = position - 1  # offsets python counting from 0
            position = position - 9  # offsets 9 bp duplication from Tn5
            tntagsList.append(position)
            i += 1
        else:
            errors += 1
    print("Reads processed " + str(i))
    print("Errors " + str(errors))

    return tntagsList
