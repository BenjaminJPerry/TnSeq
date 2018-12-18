# 2018 Benjamin J Perry - Attribution-NonCommercial-ShareAlike 4.0 International
# (CC BY-NC-SA 4.0)
# Version: 1.0.0
# Maintainer: Benjamin J Perry
# Email: benjamin.perry@postgrad.otago.ac.nz
# Status: Functional
# Citation: TBD


import argparse

def wigRefGen(fastafile, outfile, returnStatus = False):
    '''
    Reads in a single fasta sequence from a file. Identifies the position of the T in ever TA motif.
    Returns a wig files with every TA site and the count column set to 100.
    Fasta header line is used to the wig config line.

    :param fastafile: String; input fasta file.
    :param outfile: String; outfile to write reference 'TA' wig track to.
    :param returnStatus: Boolean; if True returns Ref wig file object with function call.
    :return: return list of [pos, count] of TA sites with wig header line.
    '''

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

    fastaSequence = ''
    for chunk in fastaData:
        fastaSequence = fastaSequence + str(chunk)

    # Identify the position of all TA motifs, append them to a reference list with default count of 100.
    i = 0
    TACount = 0
    TAPositions = []
    # Identify positions with TA Motifs
    while i < len(fastaSequence)-1:
        if fastaSequence[i] == 'T' and fastaSequence[i+1] == 'A':
            TACount += 1
            entry = [i, 100]
            TAPositions.append(entry)
        i += 1
        continue
    # Check for contig spanning TA site
    if fastaSequence[-1] == 'T' and fastaSequence[0] == 'A':
        TACount += 1
        entry = [fastaSequence[-1], 100]
        TAPositions.append(entry)

    # Print out the reference .wig file
    outputFile = open(outfile, 'w')
    # Print out the .wig file to outfile
    outputFile.write('variableStep chrom=' + headerLine + '\n')
    for line in TAPositions:
        outputFile.write(str(line[0]) + '\t' + str(line[1]) + '\n')
    outputFile.close()

    # Check return status option and comply
    returnList = []
    if returnStatus = True:
        returnHeader = ['variableStep chrom=' + headerLine + '\n', '']
        returnList.append(returnList)
        returnList.append(TAPositions)
        return returnList
    elif returnStatus = False:
        raise SystemExit(1)

if __name__ == '__main__':
    # Parse arguments

    # Compute Wig Reference Track
    wigRefGen(inFile, outFile)
