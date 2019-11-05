# 2018 Benjamin J Perry - Attribution-NonCommercial-ShareAlike 4.0 International
# (CC BY-NC-SA 4.0)
# Version: 2.1.0
# Maintainer: Benjamin J Perry
# Email: benjamin.perry@postgrad.otago.ac.nz
# Status: Dev
# Citation: TBD

def bedfileToTntags(bedfile, tn5=False): #TODO Re-Write to Process as Readfile vs. Holding File in Memory
    ''' ('path/to/bedfile.bed') -> tntagsList[]

    Opens, reads, and parses a bedfile given path to the file.
    Returns a list, each element the genomic coordinate of tntag start.
    This script takes into account the strand the read was mapped to.

    '''

    #Open input bedfile
    bedfileIn = open(bedfile, 'r')

    #Process input bedfile
    bedfileData = []
    for line in bedfileIn:
        items = line.rstrip('\r\n').split('\t')
        items = [item.strip() for item in items]
        bedfileData.append(items)

    #Process tntag start positions
    tntagsList = []
    i = 0
    errors = 0
    for line in bedfileData: #TODO update to check Tn5 Condition to implement the 9 bp offset
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
    '''(tntagList[], referenceWigList[][]) -> updatedWigList[][]

    Accepts a list of tn-tag positions and a reference .wig list. Sums unique tn-tag start positions and
    updates the counts for each start position in the reference wig list. Returns the updated wig list, with all
    unobserved tn-tag start positions set to 0, and all observed tn-tag start positions set to their counts.

    '''
    import collections

    updatedWig = []
    # Make Counter for position list
    positionCounter = collections.Counter(tntagList)
    # Update count for each position in referenceWig
    for siteTA in referenceWig:
        updatedEntry = []
        updatedEntry.append(siteTA[0])
        updatedEntry.append(positionCounter[siteTA[0]])
        updatedWig.append(updatedEntry)
    return updatedWig

def wigPipe(fastafile, bedfile, wigOutfile):
    ''' ("path/fastafile", "path/bedfile", "path/wigOut") -> Null; Prints output wig file for treatment.

    :param fastafile: input fasta file with a single >headerLine
    :param bedfile: tnseq data bedfile derived from read alignment
    :return: null

    Takes a fasta file, and a corresponding bedfile of tntag alignments, and return a wig file of the tntag counts.
    Reference wig files can be generated using refGenTA.wigRefGen()

    '''
    import refGenTA
    import wigScripts
    # Generate the reference wig counts list
    print('Generating reference wig position list.')
    referenceWig = refGenTA.wigRefGen(fastafile, printStatus=False)
    printHeader = referenceWig.pop(0) #refGenTA must return print header as item[0]
    # Transposon insertion counts from the bedfile
    print('\nProcessing tntag insertions bedfile.')
    tnCounts = wigScripts.bedfileToTntags(bedfile)

    # Make Update tnCount wig track #TODO Add in Stranded .wig file output files
    print('\nCreating treatment specific tntag wig track.')
    treatmentWig = wigScripts.updateWigList(tnCounts, referenceWig)

    wigPrint = open(wigOutfile, 'w')
    print('\nWriting treatment wig file to: ' + str(wigPrint))
    wigPrint.write(str(printHeader[0]) + '\n')
    printReadsCount=0
    for line in treatmentWig:
        wigPrint.write(str(line[0]) + '\t' + str(line[1]) + '\n')
        printReadsCount += line[1]
    wigPrint.close()
    print('\nTotal reads printed to .wig file: ' + str(printReadsCount))

    raise SystemExit(0)

if __name__ == '__main__':
    # Parse arguments
    import argparse
    import refGenTA
    import sys

    parser = argparse.ArgumentParser()

    if len(sys.argv[1:]) == 0:
        parser.print_usage() # for just the usage line
        parser.exit()

    parser.add_argument("-F", "--fastaFile", type = str, help="path to input .fasta file.")
    parser.add_argument("-B", "--bedFile", type = str, help="path to input .bed file.")
    parser.add_argument("-O", "--outputWigFile", type = str, help="path to output .wig file to be written.")
    args = parser.parse_args()

    inFile = args.fastaFile
    outFile = args.outputWigFile
    bedFile = args.bedFile

    # Compute Wig Reference Track
    wigPipe(inFile, bedFile, outFile)

    print("Treatment .wig file printed to: "+args.outputWigFile+"\n")

    print('2018 Benjamin J Perry - (CC BY-NC-SA 4.0)')
    print('Version: 2.1.0')
    print('Email: benjamin.perry@postgrad.otago.ac.nz')
    print('Citation: TBD\n')
