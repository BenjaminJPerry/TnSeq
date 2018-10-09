
#Creator: Benjamin Perry
#Summary: This script will randomly search through a genome and capture 16 bps downstream of 'TA'.
#         The Mariner transposon inserts on the sequence 'TA'

print "Running in silico INSeq read generator.\nThis program will capture the 16 bps upstream and downstream of a mariner insertion site selected at random.\nThe program will output the reads in .fastq format into a file titled 'in_silico_reads.txt'.\n\n\n";
print "Generating Reads...\n";

#Randomly Determine the Replicon to use
$min = 1;
$range = 100;
$loop = 0;

while ($loop <= 100){
$min = 1;
$range = 100;
$replicon = int(rand($range))+ $min;

##For Chromosome........................................................................................................
if($replicon <= 65){
print "$replicon Chromosome was selected\n";
$loop++;


$range = 2;
$strand = int(rand($range));
print " Strand $strand\n "
      #if ($strand ==0){
      #  print "- strand \n";}
      #elsif ($strand ==1){
      #   print "+ strand \n";}

}
##For pRL12........................................................................................................
elsif ($replicon > 65 && $replicon <= 76){
print "$replicon pRL12 was selected\n";
$loop++;
}
##For pRL11........................................................................................................
elsif ($replicon > 76 && $replicon <= 85){
print "$replicon pRL11 was selected\n";
$loop++;
}
##For pRL10........................................................................................................
elsif ($replicon > 85 && $replicon <= 91){
print "$replicon pRL10 was selected\n";
$loop++;
}
##For pRL9.........................................................................................................
elsif ($replicon > 91 && $replicon <= 96){
print "$replicon pRL9 was selected\n";
}
##For pRL8.........................................................................................................
elsif ($replicon > 96 && $replicon <= 98){
print "$replicon pRL8 was selected\n";
$loop++;
}
##For pRL7.........................................................................................................
elsif ($replicon <=100 && $replicon >98 ){
print "$replicon pRL7 was selected\n";
$loop++;
}
}

exit;







open (INFILE, "Rlv3841_Chr.fasta");
@rawseq = <INFILE>;
close (INFILE);

open (OUTFILE, ">>in_silico_reads.txt");
###process the array to remove carraige returns
foreach $element (@rawseq){
  chomp $element;
};

### This block could be used to record the fasta header of the file being used to generate the in silico data set.
#$header = @rawseq[0];
shift @rawseq;
#print "$header \n\n";


###joining the array into a string
$fullseq = join( '', @rawseq);

###Generating the complementary strand in 5' to 3' orientation
$revfullseq = reverse ($fullseq);
$copyseq = $revfullseq;
$copyseq =~ tr/ATCG/TAGC/;
$revcompfullseq = $copyseq;

###Defining paramters for the random number generator
$min = 1;
$range = 5057142;

### This loop will randomly select a position in the (+) strand, check if it is a 'T', if it is it will then check is the next position is an A, if it is it will then capture the 16bps up and down stream of the TA motif. It then increments the counter.
###to modify the number of in silico reads produced change the integer in the while loop logical operator.
$iteration5 =0;
while ($iteration5 < 10)### Change the integer to change the number of reads, the number of reads is equal to 2 x the number of loops.
{ ##This is were you set how many in silico insertios you want
$pos1 = int(rand($range))+ $min;
$pos2 = $pos1 + 1;
$base1 = substr $fullseq, $pos1, 1;
$base2 = substr $fullseq, $pos2, 1;

if ($base1 eq T and $base2 eq A)###checks if the position chosen is a mariner insertions site: 'TA'
{
  $capture3 = $pos1;
  $capture5 = $pos1 - 16;### To change the size of the read you capture, change all of the following integers to the size of read desired. It is currently set to capture 16 bps as this is what happens in the INSeq method.
  $read3 = substr $fullseq, $capture3, 16;
  $read5 = substr $fullseq, $capture5, 16;
  $capdown3 = $capture3+16;
  $capdown5 = $capture5+16;
  print OUTFILE "\@Read_Pos:$capture3...$capdown3\n$read3\n+\nRRRRRRRRRRRRRRRR\n";###prints the captured peice in .fastq format with a Phred+64 (Illumina1.3+) quality score of 20 (R) for each base
  print OUTFILE "\@Read_Pos:$capture5...$capdown5\n$read5\n+\nRRRRRRRRRRRRRRRR\n";###prints the captured peice in .fastq format with a Phred+64 (Illumina1.3+) quality score of 20 (R) for each base
  $iteration5++;
}

}

print "Reads generated from the (+) strand\n";
#print OUTFILE "\n\n";

###This while loop is identical to the previous loop, except it functions on the complementary (-) strand.
###Change all of the paramters in the loop to the identical values as the values in the previous loop.
$iteration3 =0;
while ($iteration3 < 10)
{
$pos1 = int(rand($range))+ $min;
$pos2 = $pos1+1;
$base1 = substr $revcompfullseq, $pos1, 1;
$base2 = substr $revcompfullseq, $pos2, 1;

if ($base1 eq T and $base2 eq A)
{
  $capture3 = $pos1;
  $capture5 = $pos1 - 16;
  $read3 = substr $revcompfullseq, $capture3, 16;
  $read5 = substr $revcompfullseq, $capture5, 16;
  $capdown3 = $capture3+16;
  $capdown5 = $capture5+16;

  print OUTFILE "\@Read_Pos:$capture3...$capdown3\n$read3\n+\nRRRRRRRRRRRRRRRR\n";
  print OUTFILE "\@Read_Pos:$capture5...$capdown5\n$read5\n+\nRRRRRRRRRRRRRRRR\n";
  $iteration3++;
}

}
print "Reads removed from the (-) strand\n";
print "Finished and exiting.\n";
exit;


#if ($base1 eq A and $base2 eq T){ ##Checling if the two positions combine to make a Mariner Insertion Site on the 3' strand

#   print "randnumb: $pos1 \nrandnumb2: $pos2\n";
#   print "Base position 1 is: $base1 \n";
#   print "Base position 2 is: $base2 \n";

#  $capture = $pos1 - 16;
#  print "$capture\n";
#  $read = substr $fullseq, $capture, 16;
#  print OUTFILE "$read\n";
#  $iteration++;
#}


