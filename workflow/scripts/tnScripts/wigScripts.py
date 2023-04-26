# 2020 Benjamin J Perry - Attribution-NonCommercial-ShareAlike 4.0 International
# (CC BY-NC-SA 4.0)
# Version: 2.3.0
# Maintainer: Benjamin J Perry
# Email: benjamin.perry@postgrad.otago.ac.nz
# Status: Dev
# Citation: TBD

def bedfileToTntags(bedfile):
    """ ('path/to/bedfile.bed') -> tntagsList[]
    :param bedfile: String; path to an input bedfile.
    :return: List[pos]; Tn insertion positions of all reads.

    Parses a bedfile and returns a list of genomic coordinate Tn insertion positions.
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
            position = position + 1
            tntagsList.append(position)
            i += 1
        elif str(line[5]) == '-':
            position = int(line[2])
            position = position - 1
            tntagsList.append(position)
            i += 1
        else:
            errors += 1
    print("Reads processed " + str(i))
    print("Errors " + str(errors))

    return tntagsList


def updateWigList(tntagList, referenceWig):
    """(tntagList[], referenceWigList[][]) -> updatedWigList[][]
    :param tntagList: list[pos]; Tn insertion positions.
    :param referenceWig: list[pos][count = 0]; Reference wig list.
    :return: list[pos][count = [insFreq]; updated reference with counts = Freq. of Pos.

    Takes a list of Tn insertions and reference wig array,
    updates reference wig array with frequency of each position in
    list of Tn insertions.
    """
    import collections

    updatedWig = []
    # Make Counter for position list
    positionCounter = collections.Counter(tntagList)
    # Update count for each position in referenceWig
    for site in referenceWig:
        updatedEntry = []
        updatedEntry.append(site[0])
        updatedEntry.append(positionCounter[site[0]])
        updatedWig.append(updatedEntry)
    return updatedWig


def wigPipe(fastafile, bedfile, wigOutfile, Tn5=False):
    """ ("path/fastafile", "path/bedfile", "path/wigOut", Tn5=False) -> Null; Prints output wig file for treatment.
    :param fastafile: String; Input reference fasta file with a single >headerLine.
    :param bedfile: String; Tnseq data bedfile derived from read alignment.
    :param wigOutfile: String; Path to write output wig file to.
    :param Tn5: Boolean; Is Tn5 Tnseq data.
    :return: null

    Takes a fasta file, and a corresponding bedfile of tntag alignments,
    returns a wig file of the tntag counts.
    Reference wig files can be generated using refGenTA.wigRefGen()
    """

    import refGenTA
    import wigScripts
    import tn5Wig

    # Generate the reference wig counts list

    if not Tn5:
        print('Generating reference wig position list.\n')
        referenceWig = refGenTA.wigRefGen(fastafile, printStatus=False)
        printHeader = referenceWig.pop(0)  # refGenTA must return print header as item[0]
        # Transposon insertion counts from the bedfile
        print('Processing tntag insertions bedfile.\n')
        tnCounts = wigScripts.bedfileToTntags(bedfile)

        # Make Update tnCount wig track
        print('Creating treatment specific tntag wig track.\n')
        treatmentWig = wigScripts.updateWigList(tnCounts, referenceWig)

        wigPrint = open(wigOutfile, 'w')
        print('Writing treatment wig file to: ' + str(wigPrint) + '\n')
        wigPrint.write(str(printHeader[0]) + '\n')
        printReadsCount = 0
        for line in treatmentWig:
            wigPrint.write(str(line[0]) + '\t' + str(line[1]) + '\n')
            printReadsCount += line[1]
        wigPrint.close()
        print('Total reads printed to .wig file: ' + str(printReadsCount) + '\n')
        raise SystemExit(0)

    if Tn5:
        print('Generating reference wig position list.\n')
        referenceWig = tn5Wig.tn5RefWig(fastafile)
        printHeader = referenceWig.pop(0)  # refWig must return print header as item[0]

        # Transposon insertion counts from the bedfile
        print('Processing tntag insertions bedfile.\n')
        tnCounts = tn5Wig.tn5BedfileToTntags(bedfile)

        # Make Update tnCount wig track
        print('Creating treatment specific tntag wig track.\n')
        treatmentWig = wigScripts.updateWigList(tnCounts, referenceWig)

        wigPrint = open(wigOutfile, 'w')
        print('Writing treatment wig file to: ' + str(wigPrint) + '\n')
        wigPrint.write(str(printHeader[0]) + '\n')
        printReadsCount = 0
        for line in treatmentWig:
            wigPrint.write(str(line[0]) + '\t' + str(line[1]) + '\n')
            printReadsCount += line[1]
        wigPrint.close()
        print('Total reads printed to .wig file: ' + str(printReadsCount) + '\n')
        raise SystemExit(0)

if __name__ == '__main__':
    # Parse arguments
    import argparse
    import sys

    parser = argparse.ArgumentParser()

    if len(sys.argv[1:]) == 0:
        parser.print_usage()  # for just the usage line
        parser.exit()

    parser.add_argument("-F", "--fastaFile", type=str, help="path to input .fasta file.")
    parser.add_argument("-B", "--bedFile", type=str, help="path to input .bed file.")
    parser.add_argument("-O", "--outputWigFile", type=str, help="path to output .wig file to be written.")
    parser.add_argument("-Tn5", "--Tn5Bedfile", action='store_true', help="Flag for Tn5 Tnseq.")
    args = parser.parse_args()

    # Parse arguments
    inFile = args.fastaFile
    outFile = args.outputWigFile
    bedFile = args.bedFile
    Tn5 = args.Tn5Bedfile

    # Compute Wig Reference Track
    wigPipe(inFile, bedFile, outFile, Tn5)

    print("Treatment .wig file printed to: " + args.outputWigFile + "\n")

    print('2020 Benjamin J Perry - (CC BY-NC-SA 4.0)')
    print('Version: 2.3.0')
    print('Email: benjamin.perry@postgrad.otago.ac.nz')
    print('Citation: TBD\n')
