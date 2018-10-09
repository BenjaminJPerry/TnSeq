## Creator: Benjamin Perry
## Summary: This script counts the number of Marinier Transposon Insertion Sites in the Rhizobium leguminosarum bv. viciae 3841 genome.
##          The Mariner Transposon will insert on either strand at the motif "TA"
##          This program search ssDNA for 'TA' and 'AT' motifs separately and combines them to determine the total number of insertion sites on both strands.


##Chromosome
$TASum = 0;

open (INFILE, "Rlv3841_Chr.fasta");
@rawseq = <INFILE>; # Reads the DNA into an array
foreach $element (@rawseq){
  chomp $element;
};

$header = @rawseq[0];
shift @rawseq;
print "$header \n\n";
$fullseq = join( '', @rawseq);
$TAcountChr =0;

while($fullseq =~ m/TA/g){$TAcountChr++};
print "The total number of mariner insertion sites in the Chromosome: $TAcountChr \n\n";
close (INFILE);

$TASum = $TASum + $TAcountChr;


## Plasmid pRL7
open (INFILE, "Rlv3841_pRL7.fasta");
@rawseq = <INFILE>;
foreach $element (@rawseq){
  chomp $element;
};

$header = @rawseq[0];
shift @rawseq;
print "$header \n\n";
$fullseq = join( '', @rawseq);
$TAcountPRL7 =0;

while($fullseq =~ m/TA/g){$TAcountPRL7++};

print "The total number of mariner insertion sites in the pRL7: $TAcountPRL7 \n\n";
close (INFILE);

$TASum = $TASum + $TAcountPRL7;

## Plasmid pRL8
open (INFILE, "Rlv3841_pRL8.fasta");
@rawseq = <INFILE>;
foreach $element (@rawseq){
  chomp $element;
};

$header = @rawseq[0];
shift @rawseq;
print "$header \n\n";
$fullseq = join( '', @rawseq);
$TAcountPRL8 =0;

while($fullseq =~ m/TA/g){$TAcountPRL8++};

print "The total number of mariner insertion sites in the pRL8: $TAcountPRL8 \n\n";
close (INFILE);

$TASum = $TASum + $TAcountPRL8;

## Plasmid pRL9
open (INFILE, "Rlv3841_pRL9.fasta");
@rawseq = <INFILE>;
foreach $element (@rawseq){
  chomp $element;
};

$header = @rawseq[0];
shift @rawseq;
print "$header \n\n";
$fullseq = join( '', @rawseq);
$TAcountPRL9 =0;

while($fullseq =~ m/TA/g){$TAcountPRL9++};

print "The total number of mariner insertion sites in the pRL9: $TAcountPRL9 \n\n";
close (INFILE);

$TASum = $TASum + $TAcountPRL9;

## Plasmid pRL10
open (INFILE, "Rlv3841_pRL10.fasta");
@rawseq = <INFILE>;
foreach $element (@rawseq){
  chomp $element;
};

$header = @rawseq[0];
shift @rawseq;
print "$header \n\n";
$fullseq = join( '', @rawseq);
$TAcountPRL10 =0;

while($fullseq =~ m/TA/g){$TAcountPRL10++};

print "The total number of mariner insertion sites in the pRL10: $TAcountPRL10 \n\n";
close (INFILE);

$TASum = $TASum + $TAcountPRL10;

## Plasmid pRL11
open (INFILE, "Rlv3841_pRL11.fasta");
@rawseq = <INFILE>;
foreach $element (@rawseq){
  chomp $element;
};

$header = @rawseq[0];
shift @rawseq;
print "$header \n\n";
$fullseq = join( '', @rawseq);
$TAcountPRL11 =0;

while($fullseq =~ m/TA/g){$TAcountPRL11++};

print "The total number of mariner insertion sites in the pRL11: $TAcountPRL11 \n\n";
close (INFILE);

$TASum = $TASum + $TAcountPRL11;

## Plasmid pRL12
open (INFILE, "Rlv3841_pRL12.fasta");
@rawseq = <INFILE>;
foreach $element (@rawseq){
  chomp $element;
};

$header = @rawseq[0];
shift @rawseq;
print "$header \n\n";
$fullseq = join( '', @rawseq);
$TAcountPRL12 =0;

while($fullseq =~ m/TA/g){$TAcountPRL12++};

print "The total number of mariner insertion sites in the pRL12: $TAcountPRL12 \n\n";
close (INFILE);

$TASum = $TASum + $TAcountPRL12;
print "The total number of mariner transposon insertions sites is: $TASum.";
exit;
