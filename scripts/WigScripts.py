# 2018 Benjamin J Perry - Attribution-NonCommercial-ShareAlike 4.0 International
# (CC BY-NC-SA 4.0)
# Version: 0.0.1
# Maintainer: Benjamin J Perry
# Email: benjamin.perry@postgrad.otago.ac.nz
# Status: Dev
# Citation: TBD

def bedfileToTntags(bedfile): # TODO Validate Function if Working
    ''' ('path/to/bedfile.bed') -> tntagsList[]

    Opens, reads, and parses a bedfile given path to the file.
    Returns a list, each element the genomic coordinate of tntag 'T' start.
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
    print "Reads processed " + str(i)
    print "Errors " + str(errors)

    return tntagsList

def updateWigList(tntagList, referenceWig): # TODO Validate Function if Working
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
        updatedEntry.append(int(siteTA[0]))
        updatedEntry.append(int(positionCounter[int(siteTA[0])]))
        updatedWig.append(updatedEntry)
