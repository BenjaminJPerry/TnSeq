### Count the number of genes in specific functional classes
###
###

### Read in the .csv file
#User defines .csv input file
file1=raw_input('Designate the input locus_state.csv file: ')
print "File:", file1
csv_in = open(file1, 'r')

csv_data=[] ### 2D Array of the TnSeq Data ###
for line in csv_in:
    items = line.rstrip('\r\n').split(',')
    items = [item.strip() for item in items]
    csv_data.append(items)


### Reads in the Riley_Keys .csv file
#User defines .csv riley_funct input file
file2=raw_input('Designate the input riley_key.csv: ')
print "File:", file2
riley_in = open(file2, 'r')
#Read in locus_funct data into a 2D list

riley_key=[] ### 2D Array of the Riley Functional Class Keys ###

for line in riley_in:
    items = line.rstrip('\r\n').split(',')
    items = [item.strip() for item in items]
    riley_key.append(items)

### Makes a dictionary of the riley_key data. Riley Class number as key, string describing the class returned.
riley={}
key={}

for i in riley_key:
    key[i[0]]=i[1]

core_genome={}
core_total=0
for i in riley_key:
    core_genome[i[0]] = 0
    
vmm_unique={}
vmm_total=0
for i in riley_key:
    vmm_unique[i[0]] = 0
    
ty_unique={}
ty_total=0
for i in riley_key:
    ty_unique[i[0]] = 0

    
### Start the counting process for TY and VMM
for line in csv_data:
    if (line[8] == 'ES' or line[8] == 'GD') and (line[16] == 'ES' or line[16] == 'GD'):
        core_genome[line[17]]=core_genome[line[17]]+1
        core_total=core_total+1
        
for line in csv_data:
    if (line[8] == 'ES' or line[8] == 'GD') and (line[16] == 'NE' or line[16] == 'GA'):
        vmm_unique[line[17]]=vmm_unique[line[17]]+1
        vmm_total=vmm_total+1

for line in csv_data:
    if (line[8] == 'NE' or line[8] == 'GA') and (line[16] == 'ES' or line[16] == 'GD'):
        ty_unique[line[17]]=ty_unique[line[17]]+1
        ty_total=ty_total+1

### Print output files with the summary of functional groupings
file3 = '/home/benji/Desktop/core_genome.csv'
file4 = '/home/benji/Desktop/vmm_unique.csv'
file5 = '/home/benji/Desktop/ty_unique.csv'

outfile1 = open(file3, 'w')
outfile2 = open(file4, 'w')
outfile3 = open(file5, 'w')

core_out=[]
vmm_out=[]
ty_out=[]

for i in core_genome.items():
    funct=key[i[0]]
    j=list(i)
    j.append(funct)
    line = ",".join(str(element) for element in j)
    core_out.append(line)

for i in vmm_unique.items():
    funct=key[i[0]]
    j=list(i)
    j.append(funct)
    line = ",".join(str(element) for element in j)
    vmm_out.append(line)

for i in ty_unique.items():
    funct=key[i[0]]
    j=list(i)
    j.append(funct)
    line = ",".join(str(element) for element in j)
    ty_out.append(line)

for i in core_out:
    outfile1.write(str(i))
    outfile1.write('\n')
outfile1.close()

for i in vmm_out:
    outfile2.write(str(i))
    outfile2.write('\n')
outfile2.close()

for i in ty_out:
    outfile3.write(str(i))
    outfile3.write('\n')
outfile3.close()
