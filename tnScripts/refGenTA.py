# 2019 Benjamin J Perry - Attribution-NonCommercial-ShareAlike 4.0 International
# (CC BY-NC-SA 4.0)
# Version: 2.1.0
# Maintainer: Benjamin J Perry
# Email: benjamin.perry@postgrad.otago.ac.nz
# Status: Functional
# Citation: TBD

def wigRefGen(fastafile, outfile='', printStatus=False):
    '''(reference.fasta, output.wig, printStatus=Boolean) -> refWig[][] OR Null; print to output.wig if True
    Reads in a single fasta sequence from a file. Identifies the position of the T in ever TA motif.
    Returns a wig files with every TA site and the count column set to 100.
    Fasta header line is used to the wig config line.

    :param fastafile: String; input fasta file.
    :param outfile: String; outfile to write reference 'TA' wig track to.
    :param printStatus: Boolean; if True prints Ref wig file object with function call.
    :return: return list of [pos, count] of TA sites with wig header line.
    '''
    #TODO Implement Checks for multi-fasta inputs

    print('\nComputing Reference TA .wig file for: ' + fastafile)
    inputFile = open(fastafile, 'r')
    fastaData = inputFile.read()

    # Check fastafile for a single header line.
    if fastaData.count('>') > 1:
        print('Error: More than one fasta header line in ' + fastafile)
        raise SystemExit(1)

    # Process the fasta file into a header and a sequence variables
    fastaData = fastaData.splitlines()
    headerLine = fastaData.pop(0)
    headerLine = headerLine.split()[0]
    headerLine = headerLine[1:]

    fastaSequence = ''
    for chunk in fastaData:
        fastaSequence = fastaSequence + str(chunk)

    # Identify the position of all TA motifs, append them to a reference list with default count of 100.
    i = 0
    TACount = 0
    TAPositions = []
    printHeader = ['variableStep chrom=' + str(headerLine)]
    TAPositions.append(printHeader)

    # Identify positions with TA Motifs
    while i < len(fastaSequence) - 1:
        if fastaSequence[i] == 'T' and fastaSequence[i + 1] == 'A':
            TACount += 1
            entry = [i + 1, 100]
            TAPositions.append(entry)
        i += 1
        continue

    # Check for contig spanning TA site
    if fastaSequence[-1] == 'T' and fastaSequence[0] == 'A':
        TACount += 1
        entry = [fastaSequence[-1], 100]
        TAPositions.append(entry)

    print('Total TA motifs found: ' + str(TACount))

    # Check print status option and comply
    if not printStatus:
        return TAPositions

    if printStatus:
        print('Reference TA .wig file computed.')
        print('Printing reference TA .wig file.')
        outputFile = open(outfile, 'w')
        # Print out the .wig file to outfile
        printHeader = TAPositions.pop(0)
        printHeader = str(printHeader[0])
        outputFile.write(printHeader + '\n')
        for line in TAPositions:
            outputFile.write(
                str(line[0]) + '\t' + str(line[1]) + '\n')
        outputFile.close()


if __name__ == '__main__':
    # Parse arguments
    import argparse
    import sys

    parser = argparse.ArgumentParser()

    if len(sys.argv[1:]) == 0:
        parser.print_usage()
        parser.exit()

    parser.add_argument("-F", "--fastaFile", type=str, help="path to input .fasta file.")
    parser.add_argument("-O", "--outputWigFile", type=str, help="path to output .wig file to be written.")
    args = parser.parse_args()

    inFile = args.fastaFile
    outFile = args.outputWigFile

    # Compute Wig Reference Track
    wigRefGen(inFile, outfile=outFile, printStatus=True)

    print("Reference .wig file printed to: " + args.outputWigFile + "\n")

    print('2018 Benjamin J Perry - (CC BY-NC-SA 4.0)')
    print('Version: 2.1.0')
    print('Email: benjamin.perry@postgrad.otago.ac.nz')
    print('Citation: TBD\n')
