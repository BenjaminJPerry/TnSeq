#Creator: Benjamin Perry
#Summary: This Script will read a .gff file and the coresponding .fasta file and count the number of TA mariner
#	  Insertion Sites in each ORF. It will then out put the locus and TA insertion number in a tab delimited format.
#	  It will also output each insertion in each ORF and what percentage of the ORF the TA insertion is located in into a seperate file


## Opening the Fasta Sequence File, Reading it into an array
open (INFILE, "H_INF.fasta"); ### <<<<<<<<<<<<<<<<<<<<<<<<< EDIT .fasta SEQUENCE INPUT FILE HERE
@rawseq = <INFILE>;
close (INFILE);
open (OUTFILE1, ">>HINF_ORF_TA_FREQ.txt");### <<<<<<<<<<<<<<<< EDIT OUTPUT FOR ORF TA FREQUENCY DATA FILE HERE
print OUTFILE1 "Locus	Start	End	Strand	Length	TA_Ins\n";

open (OUTFILE2, ">>HINF_ORF_TA_REL_POS.txt");### <<<<<<<<<<<<<<<<< EDIT OUTPUT FOR ORF TA REL POS DATA FILE HERE
print OUTFILE2 "Locus	Ins_Pos\n";

## Removing Carraige Returns
foreach $section (@rawseq)
{
	chomp $section;
};
shift @rawseq;

##joining the array into a string
$fullseq = join( '', @rawseq);

## Opening the .gff file and processing each feature
open(INFILE,"H_INF.gff"); ### <<<<<<<<<<<<<<<<<<<<<<<<<<<< EDIT .gff FEATURE INPUT FILE HERE
while ( <INFILE> )
{
chomp;
@feature = split /\s/;

## SubString Variables
$condition = @feature[2];
$start = @feature[3];
$startSS = $start-1;
$end   = @feature[4];
$length= $end-$start;
$strand= @feature[6];
$TAInst= 0;

## Break apart the feature info to collect the locus tag
$FeatInfo = @feature[8];
@meta = split(/;/, $FeatInfo);
@locustag = split( /=/, @meta[0]);
$locus = @locustag[1];


## Capture the coding sequence of the .gff feature
$ORF = substr $fullseq, $startSS, $length;

##Filter the info
if($condition eq "CDS")
{

## Check to see the strand the ORF is located on; Makes $ORF the reverse compliment of the + strand
if($strand eq "-")
{
$revORF = reverse($ORF);
$revORF =~ tr/ATCG/TAGC/;
$ORF    = $revORF;
$i = $start;
$start = $end;
$end = $i+1;
} 

### Now we count the TA Insertions in the ORF
$ORFCount = $ORF;
while($ORFCount =~ m/TA/g){$TAInst++};


########### Output into a tab delimited .txt file
print OUTFILE1 "$locus	$start	$end	$strand	$length	$TAInst\n";	

### Now we calculate what percentage of the sequence, relative to the start codon, the TA is in
$TAPerc = $ORF;
$TA = "TA";
$offset = 0;
$pos = index($TAPerc, $TA, $offset);
while($pos != -1)
{
$PosPercent = ($pos/$length)*100;
print OUTFILE2 "$locus	";
printf OUTFILE2 ("%.1f ", $PosPercent);
print OUTFILE2 "\n";
$offset = $pos+1;
$pos = index($TAPerc, $TA, $offset);
};
}	

###################################################
}
exit;
