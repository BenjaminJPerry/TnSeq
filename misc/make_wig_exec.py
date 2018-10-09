#!/usr/bin/env python
import collections
import os
import getopt
import sys

passed=str(sys.argv)
wkdir=os.getcwd()
file1=wkdir+'/'+str(sys.argv[1])
print "File:", file1
bedfile = open(file1, 'r')

# Read bedfile into list
data=[]

for line in bedfile:
    items = line.rstrip('\r\n').split('\t')
    items = [item.strip() for item in items]
    data.append(items)
#Parse the reads in the bedfile into replicon specific lists
datachr=[]
data12=[]
data11=[]
data10=[]
data09=[]
data08=[]
data07=[]
bedfile_errors=[]
for line in data:
        if str(line[0]) == 'AM236080': #Chrom
            datachr.append(line)            
        elif str(line[0]) == 'AM236086': #pRL12
            data12.append(line)
        elif str(line[0]) == 'AM236085': #pRL11
            data11.append(line)
        elif str(line[0]) == 'AM236084': #pRL10
            data10.append(line)
        elif str(line[0]) == 'AM236083': #pRL09
            data09.append(line)
        elif str(line[0]) == 'AM236082': #pRL08
            data08.append(line)
        elif str(line[0]) == 'AM236081': #pRL07
            data07.append(line)
        else:
            bedfile_errors.append(line)
# Append each replicon data file to a list
allreads=[]
allreads.append(datachr)
allreads.append(data12)
allreads.append(data11)
allreads.append(data10)
allreads.append(data09)
allreads.append(data08)
allreads.append(data07)

# Loop over allreads, find 'T' positons, append positions to a list allTs
allTs=[]
errors=0
for bedfile in allreads:
    tempTs=[]
    for line in bedfile:
        if str(line[5]) == '+':
            position = int(line[1])
            position = position+1
            tempTs.append(position)
        elif str(line[5]) == '-':
            position = int(line[2])
            position = position-1
            tempTs.append(position)
        else: errors=errors+1
    allTs.append(tempTs)

# Read in and prepare reference wig files
chrwig='/home/benji/Dropbox/Rlv3841_Reference_Files/AM236080_Reference.wig'
wigchrom = open(chrwig, 'r')
chromref=[]
for line in wigchrom:
    items = line.rstrip('\r\n').split(' ')
    items = [item.strip() for item in items]
    chromref.append(items)
    
chromref=chromref[2:]
for line in chromref:
    line[1]=0
    
wigfile12='/home/benji/Dropbox/Rlv3841_Reference_Files/AM236086_Reference.wig'
wig12 = open(wigfile12, 'r')
pRL12ref=[]
for line in wig12:
    items = line.rstrip('\r\n').split(' ')
    items = [item.strip() for item in items]
    pRL12ref.append(items)
    
pRL12ref=pRL12ref[2:]
for line in pRL12ref:
    line[1]=0
       
wigfile11='/home/benji/Dropbox/Rlv3841_Reference_Files/AM236085_Reference.wig'
wig11 = open(wigfile11, 'r')
pRL11ref=[]
for line in wig11:
    items = line.rstrip('\r\n').split(' ')
    items = [item.strip() for item in items]
    pRL11ref.append(items)
    
pRL11ref=pRL11ref[2:]
for line in pRL11ref:
    line[1]=0
        
wigfile10='/home/benji/Dropbox/Rlv3841_Reference_Files/AM236084_Reference.wig'
wig10 = open(wigfile10, 'r')
pRL10ref=[]
for line in wig10:
    items = line.rstrip('\r\n').split(' ')
    items = [item.strip() for item in items]
    pRL10ref.append(items)
    
pRL10ref=pRL10ref[2:]
for line in pRL10ref:
    line[1]=0
        
wigfile9='/home/benji/Dropbox/Rlv3841_Reference_Files/AM236083_Reference.wig'
wig09 = open(wigfile9, 'r')
pRL09ref=[]
for line in wig09:
    items = line.rstrip('\r\n').split(' ')
    items = [item.strip() for item in items]
    pRL09ref.append(items)
    
pRL09ref=pRL09ref[2:]
for line in pRL09ref:
    line[1]=0
        
wigfile8='/home/benji/Dropbox/Rlv3841_Reference_Files/AM236082_Reference.wig'
wig08 = open(wigfile8, 'r')
pRL08ref=[]
for line in wig08:
    items = line.rstrip('\r\n').split(' ')
    items = [item.strip() for item in items]
    pRL08ref.append(items)
    
pRL08ref=pRL08ref[2:]
for line in pRL08ref:
    line[1]=0
        
wigfile7='/home/benji/Dropbox/Rlv3841_Reference_Files/AM236081_Reference.wig'
wig07 = open(wigfile7, 'r')
pRL07ref=[]
for line in wig07:
    items = line.rstrip('\r\n').split(' ')
    items = [item.strip() for item in items]
    pRL07ref.append(items)
    
pRL07ref=pRL07ref[2:]
for line in pRL07ref:
    line[1]=0
    
# Append the wig reference lists together
wigreferences=[]
wigreferences.append(chromref)
wigreferences.append(pRL12ref)
wigreferences.append(pRL11ref)
wigreferences.append(pRL10ref)
wigreferences.append(pRL09ref)
wigreferences.append(pRL08ref)
wigreferences.append(pRL07ref)

# Loop through each replicon T position file and wig reference file

i=0
updatedwigs=[]
for replicon in allTs:
    locallist=[]
    # Make Counter for position list
    positioncount=collections.Counter(replicon)
    for refT in wigreferences[i]:
        #refpos=refT[0]
        tempcount=0
        tempcount=positioncount[int(refT[0])]
        templist=[]
        templist.append(int(refT[0]))
        templist.append(tempcount)
        locallist.append(templist)
    updatedwigs.append(locallist)
    i=i+1
    
# Print out updated INSeq wig files
replicon_names=['AM236080','AM236086','AM236085','AM236084','AM236083','AM236082','AM236081']
namedex=0
for replicon in updatedwigs:
    wkdir=os.getcwd()
    output_wig = open(wkdir+'/'+replicon_names[namedex]+'.wig', 'w')
    output_wig.write('variableStep\n')
    for line in replicon:
        output_wig.write(str(line[0]) + ' ' + str(line[1]) + '\n')
    output_wig.close()
    namedex=namedex+1
    print "Data printed to:" + str(output_wig)

exit()



    
