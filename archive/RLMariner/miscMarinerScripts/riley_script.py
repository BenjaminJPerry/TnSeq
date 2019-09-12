####Reads in a locus_state file, a locus_function annotation file, and outputs
####which grouping each locus belongs to, based on the groupings in the 3841 genome
####paper.

### Read in the Locus_State .csv file
#User defines .csv Locus_State input file
file1=raw_input('Designate the input locus_state.csv file: ')
print "File:", file1
locus_state_in = open(file1, 'r')
#Read in locus_state data into a 2D list

locus_state=[] ### 2D Array of the TnSeq Data ###

for line in locus_state_in:
    items = line.rstrip('\r\n').split(',')
    items = [item.strip() for item in items]
    locus_state.append(items)

### Read in the Locus_Func .csv file
#User defines .csv Locus_Func input file
file2=raw_input('Designate the input locus_func.csv file: ')
print "File:", file2
locus_func_in = open(file2, 'r')
#Read in locus_state data into a 2D list

locus_func=[] ### 2D Array of the TnSeq Data ###

for line in locus_func_in:
    items = line.rstrip('\r\n').split(',')
    items = [item.strip() for item in items]
    locus_func.append(items)

### Reads in the Riley_Keys .csv file
#User defines .csv riley_funct input file
file3=raw_input('Designate the input riley_key.csv: ')
print "File:", file3
riley_in = open(file3, 'r')
#Read in locus_funct data into a 2D list

riley_key=[] ### 2D Array of the Riley Functional Class Keys ###

for line in riley_in:
    items = line.rstrip('\r\n').split(',')
    items = [item.strip() for item in items]
    riley_key.append(items)

#########################################################################################
### Data has been read into locus_state[], locus_func[], and riley_key[] ################
#########################################################################################

### Makes a dictionary of the riley_key data. Riley Class number as key, string describing the class returned.
riley={}
key_check=''
for i in riley_key:
    riley[i[0]] = i[1]
    key_check=key_check+i[0]

### Make a dictionary of the locus_funct data. Locus being the key, function number returned.
func={}
for i in locus_func:
    func[i[0]] = i[1]

### Check Condition Setup for Querying the Dictionary
locus_string=''
for i in locus_func:
    locus_string=locus_string+i[0]

### Define a global list to append the appended entries too
compiled_data=[]
for i in locus_state:
    if i[0] in locus_string:
        funct_key=func[i[0]]
        i.append(funct_key)
        compiled_data.append(i)

### Now find the functional category for the key
for i in compiled_data:
    if i[2] in key_check:
        function=riley[i[2]]
        i.append(function)

##### Merge compiled_data[2] into a single integer for comparison
##for i in compiled_data:
##    key = i[2]
##    split = key.split('.')
##    merged=''
##    for j in split:
##        merged=merged+j
##    if len(merged) > 3:
##        merged=merged[0:3]
##    finished_key=int(merged)
##    i.append(finished_key)

##### Now we check each gene for the class it belongs to
##for i in compiled_data:
##    if i[4] < 1:
##        i.append('No Known Function')
##        i.append('1')
##    elif  1 <= i[4] <= 2:
##        i.append('Conserved Hypothetical')
##        i.append('2')
##    elif 100 <= i[4] <= 131:
##        i.append('Cell Processes')
##        i.append('3')
##    elif 140 <= i[4] <= 144:
##        i.append('Protection Responses')
##        i.append('4')
##    elif 150 <= i[4] <= 155:
##        i.append('Transport/Binding Proteins')
##        i.append('5')
##    elif 160 <= i[4] <= 164:
##        i.append('Adaptation')
##        i.append('6')    
##    elif i[4] == 171:
##        i.append('Cell Division')
##        i.append('7')
##    elif 200 <= i[4] <= 214:
##        i.append('Macromolecule Metabolism')
##        i.append('8')
##    elif 220 <= i[4] <= 229:
##        i.append('Macromolecule Synthesis')
##        i.append('9')
##    elif 300 <= i[4] <= 319:
##        i.append('Metabolism of Small Molecules')
##        i.append('10')
##    elif 320 <= i[4] <= 329:
##        i.append('Biosynthesis of Co-factors and Carriers')
##        i.append('11')
##    elif 330 <= i[4] <= 339:
##        i.append('Central Intermediate Metabolism')
##        i.append('12')
##    elif 340 <= i[4] <= 345:
##        i.append('Degradation of Small Molecules')
##        i.append('13')      
##    elif 350 <= i[4] <= 359:
##        i.append('Energy and Carbon Metabolism')
##        i.append('14')
##    elif 360 <= i[4] <= 361:
##        i.append('Fatty Acid Biosynthesis')
##        i.append('15')
##    elif 370 <= i[4] <= 372:
##        i.append('Nucleotide Biosynthesis')
##        i.append('16')
##    elif 400 <= i[4] <= 419:
##        i.append('Cell Envelope')
##        i.append('17')
##    elif 420 <= i[4] <= 423:
##        i.append('Ribosome Constituents')
##        i.append('18')
##    elif 500 <= i[4] <= 515:
##        i.append('Foreign DNA')
##        i.append('19')
##    elif 600 <= i[4] <= 650:
##        i.append('Regulation')
##        i.append('20')

###Print to an output file the updated .csv file
file5 = '/home/benji/Desktop/Riley_Calss_Output.out'
outfile = open(file5, 'w')

print_out=[]
for i in compiled_data:
        line = ",".join(str(element) for element in i)
        print_out.append(line)

for i in print_out:
    outfile.write(str(i))
    outfile.write('\n')
outfile.close()
