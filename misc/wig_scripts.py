def wig_pipe(infile, outfile):
    ''' (infile, outfile)

    Reads in a .bedfile formatted alignment file and returns a .wig formatted data file.

    '''
    import wig_scripts
    bedfile = wig_scripts.read_tab_data(infile)
    read_pos = wig_scripts.bedfile_to_read_pos(bedfile)
    data = wig_scripts.count_reads(read_pos)
    wig_scripts.wig_out(outfile, data)
    
    
def read_tab_data(infile):
    ''' (file) -> list[][]

    Will read in a data file that is opended and
    tab delimited, and return a 2D list with each element list[i] representing a line
    and each element list[i][j] representing the tab-delimited data in that
    line. The script will remove any whitespace characters attached to the data.

    >>>data = read_tab_data(infile)
    print data

    [['value01', value02, value03'], [value11, value12, value13], ... [valuenj, valuenj+1, valuenj+2]]
    '''

    data=[]
    for line in infile:
        items = line.rstrip('\r\n').split('\t')
        items = [item.strip() for item in items]
        data.append(items)

    return data

def read_space_data(infile):
    ''' (file) -> list[][]

    Will read in a data file that is opended and
    space delimited, and return a 2D list with each element list[i] representing a line
    and each element list[i][j] representing the tab-delimited data in that
    line. The script will remove any whitespace characters attached to the data.

    >>>data = read__space_data(infile)
    print data

    [['value01', value02, value03'], [value11, value12, value13], ... [valuenj, valuenj+1, valuenj+2]]
    '''

    data=[]
    for line in infile:
        items = line.rstrip('\r\n').split(' ')
        items = [item.strip() for item in items]
        data.append(items)

    return data

def bedfile_to_read_pos(bedfile):
    ''' (data[][]) -> read_positions[]

    Reads a 2D list form of a bedfile, and returns an array with each element containing the position of the 'T'
    in the mariner 'TA' motif. This script takes into account the strand the read was mapped to.
        
    '''
    reads = []
    i=0
    errors=0
    for line in bedfile:
        if str(line[5]) == '+':
            position = int(line[1])
            position = position+1
            reads.append(position)
            i=i+1
        elif str(line[5]) == '-':
            position = int(line[2])
            position = position-1
            reads.append(position)
            i=i+1
        else: errors=errors+1
    print "Reads processed " + str(i)
    print "Errors " + str(errors)
    
    return reads
    
def count_reads(read_pos):
    ''' (list[]) - > wig[][]

    Reads in a list of position data from bedFile_to_read_pos(). Counts the number of
    repeats that occur at a position (i.e. representing the read depth at that position)
    and outputs a 2D list with list[position][read_depth].

    '''
    repeats = []
    output = []
    repeats =[]
    processed = 0
    
    for pos in read_pos:
        if pos in repeats:
            processed = processed+1
        else:
            repeats.append(pos)
            read_number = read_pos.count(pos)
            items =[]
            items = [pos, read_number]
            output.append(items)
            processed = processed+1

    return output

def wig_out(outfile, data):
    ''' (data[][] -> .wig)

    Takes an outfile stream and a wig[][] data list as arguements and prints the data into
    the designated file in .WIG format.

    '''
    
    outfile.write('variableStep chrom= \n')
    for line in data:
        outfile.write(str(line[0]) + '	' + str(line[1]) + '\n')

    outfile.close()

    print "Data printed to:" + str(outfile)
    
def parse_bed(infile, directory):
    ''' (infile, directory) -> Print to parsed Output Files in directory

    Reads a bedfile, parses the data based on the chrom identifier.
    Prints reads to output files with the file name of chrom identifier "NC_0083XX.1.bed" onto
    the desktop.

    NOTE: This program is Rhizobium leguminosarum bv. viciae 3841 specific
    Possible chrom identifiers are:
    
    gi|116248676|ref|NC_008378.1|
    gi|116249460|ref|NC_008379.1|
    gi|116249766|ref|NC_008380.1|
    gi|116254467|ref|NC_008381.1|
    gi|116254910|ref|NC_008382.1|
    gi|116255067|ref|NC_008383.1|
    gi|116255200|ref|NC_008384.1|

    '''

    import wig_scripts
    data = wig_scripts.read_tab_data(infile)
    dir_out = directory
    
    f78 = open(dir_out+'/NC_008378.1.bed', 'w')
    f79 = open(dir_out+'/NC_008379.1.bed', 'w')
    f80 = open(dir_out+'/NC_008380.1.bed', 'w')
    f81 = open(dir_out+'/NC_008381.1.bed', 'w')
    f82 = open(dir_out+'/NC_008382.1.bed', 'w')
    f83 = open(dir_out+'/NC_008383.1.bed', 'w')
    f84 = open(dir_out+'/NC_008384.1.bed', 'w')
    err = open(dir_out+'/Errors.txt', 'w')

    for line in data:
        if str(line[0]) == 'gi|116248676|ref|NC_008378.1|':
            for i in line:
                f78.write(str(i)+'\t')
            f78.write('\n')
        elif str(line[0]) == 'gi|116249460|ref|NC_008379.1|':
            for i in line:
                f79.write(str(i)+'\t')
            f79.write('\n')
        elif str(line[0]) == 'gi|116249766|ref|NC_008380.1|':
            for i in line:
                f80.write(str(i)+'\t')
            f80.write('\n')
        elif str(line[0]) == 'gi|116254467|ref|NC_008381.1|':
            for i in line:
                f81.write(str(i)+'\t')
            f81.write('\n')
        elif str(line[0]) == 'gi|116254910|ref|NC_008382.1|':
            for i in line:
                f82.write(str(i)+'\t')
            f82.write('\n')
        elif str(line[0]) == 'gi|116255067|ref|NC_008383.1|':
            for i in line:
                f83.write(str(i)+'\t')
            f83.write('\n')
        elif str(line[0]) == 'gi|116255200|ref|NC_008384.1|':
            for i in line:
                f84.write(str(i)+'\t')
            f84.write('\n')
        else:
            for i in line:
                err.write(str(i)+'\t')
            err.write('\n')
            
    f78.close()
    f79.close()
    f80.close()
    f81.close()
    f82.close()
    f83.close()
    f84.close()
    err.close()

def process_bedfile_for_HMM(infile1, infile2, outfile1):
    ''' (total_TA_pos.wig, TnSeq_read_pos.bedfile, tnseq_data.wig) -> write to tnseq_data.wig

    This program reads in a .wig file of the 'TA' positions in the Rlv3841 replicon in question,
    a bedfile of all reads from a TnSeq alignment, and combines the positon data from the bedfile
    with the .wig list. It then prints out a file with the TnSeq data in a format useable by
    the 
    
    '''
    import wig_scripts
    output = []
    processed = 0
    total_ta_pos = []
    total_tnseq_pos = []

    TA_in = wig_scripts.read_space_data(infile1)
    TnSeq_in = wig_scripts.read_tab_data(infile2)
    read_counts = wig_scripts.bedfile_to_read_pos(TnSeq_in)

    header = TA_in.pop(0)
    ### Takes the position of every 'TA' and pushes it into a list in int format
    for pos in TA_in:
        ta = int(pos[0])
        total_ta_pos.append(ta)

    for read in read_counts:
        value = int(read)
        total_tnseq_pos.append(value)
      
    for ta_site in total_ta_pos:
        read_number = total_tnseq_pos.count(ta_site)
        items =[]
        items = [ta_site, read_number]
        output.append(items)
        processed = processed+1

    outfile1.write(header[0] + ' ' + header[1] + '\n')

    for result in output:
        position = str(result[0])
        count = str(result[1])
        outfile1.write(position)
        outfile1.write(' ')
        outfile1.write(count)
        outfile1.write('\n')

    outfile1.close()
    infile1.close()
    infile2.close()
    
    
