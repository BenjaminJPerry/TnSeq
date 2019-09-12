import sys
import wig_scripts
import os

#User defines .bed input file
file1=raw_input('Designate the input .bed formatted file: ')
print "File:", file1
raw_bed = open(file1, 'r')

#User defines a working directory
dir_out=raw_input('Designate the working directory: ')
print "Selected: ", dir_out

os.makedirs(dir_out+'/3_conversion')
bed_dir = dir_out+'/3_conversion'

#Parse the .bed file into a working directory (modify where in the module)
wig_scripts.parse_bed(raw_bed, bed_dir)

#loop through each replicon:
#1) wig_out -> a replicon unique collection point
os.makedirs(dir_out+'/4_hmm')
hmm_dir = dir_out+'/4_hmm'

loop=78
while loop < 85:
    bed_file = bed_dir+'/NC_0083'+str(loop)+'.1.bed'
    wig_out = bed_dir+'/NC_0083'+str(loop)+'.wig'
    bed_in = open(bed_file, 'r')
    wig_o = open(wig_out, 'w')

    wig_scripts.wig_pipe(bed_in, wig_o)
    
#2) prepare .wig for HMM analysis
    hmm_out = hmm_dir+'/NC_0083'+str(loop)+'.hmm.wig'
    ref_wig = '/usr/lib/tnseq/wig_ref/ref_NC_0083'+str(loop)+'.wig'

    bed_in = open(bed_file, 'r')
    ref_in = open(ref_wig, 'r')
    hmm_wig = open(hmm_out, 'w')

    wig_scripts.process_bedfile_for_HMM(ref_in, bed_in, hmm_wig)

    loop=loop+1

sys.exit()
