def tntag_clean(infile, outfile):
        '''(.fastq -> .fasta)

        Reads in a .fastq file of tntags with adapters trimmed using cutadapt,
        checks them for 'TA' motif at the begining and that length is less then
        16bp and greater then 14bp. If greater than 16bp, the tag is trimmed to
        16bp in length. Likely due to an insertion in the tntag or to do with
        the adapter ligation in library preparation.

        '''
        raw_reads=[]
        
        #Read in the raw fastq file, each line of the file stored as a list element

        with infile as f:
                raw_reads = f.read().splitlines()

        print 'Reading raw reads...\n'
        reads = (len(raw_reads)/4)
        print reads, 'Reads read...\n'

        #Check each TnTag for two conditions and push the tag into tntags[].
        
        tntag=''
        trim_tags=[]
        trimmed=1
        count=1
        QC_failed=0
        tntags=[]
        
        for tag in raw_reads:
                tntag = tag
                if len(tntag) > 13 and len(tntag) < 17:
                        if tntag[0] == 'T' and tntag[1] == 'A':
                                tntags.append(tntag)
                                count=count+1
                        else:
                                QC_failed = QC_failed+1
                elif len(tntag) > 16:
                        if tntag[0] == 'T' and tntag[1] == 'A':
                                tntag = tntag[0:15]
                                tntags.append(tntag)
                                trimmed=trimmed+1
                        else:
                                QC_failed = QC_failed+1
         
        #Now we push the tntags into an output file.
        read_num = 1
        for tntag in tntags:
                outfile.write('>read'+str(read_num)+'\n')
                outfile.write(tntag+'\n')
                int(read_num)
                read_num = read_num + 1

        print count, 'reads had adapter sequence...\n'
        print trimmed, 'reads were trimmed...\n'        
        print read_num, 'processed tntags printed to outputfile', outfile

