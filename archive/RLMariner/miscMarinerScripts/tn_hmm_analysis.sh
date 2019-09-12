#!/bin/bash

cd "/home/benji/"

echo `python tn_hmm.py`

printf "Enter the working directory:"
read WORKDIR

echo`mkdir $WORKDIR/5_results`

RESULTDIR=$WORKDIR/5_results
HMMDIR=$WORKDIR/4_hmm
#echo "$RESULTDIR"
#echo "$HMMDIR"
PREFIX="/NC_0083"
ISUFFIX=".hmm.wig"
OSUFFIX=".out"
RSUFFIX=".result.txt"
GFFPRE="/usr/lib/tnseq/gff/NC_0083"
GFFSUF=".1_v2.gff3"

I=78

while [ $I -lt 85 ]; do
   echo "Input file: $HMMDIR$PREFIX$I$ISUFFIX"
   echo "Intermediate file: $HMMDIR$PREFIX$I$OSUFFIX"
   echo "Output file: $HMMDIR$PREFIX$I$RSUFFIX"
    INFILE=$HMMDIR$PREFIX$I$ISUFFIX
    INTFILE=$HMMDIR$PREFIX$I$OSUFFIX
    OUTFILE=$RESULTDIR$PREFIX$I$RSUFFIX
    GFF=$GFFPRE$I$GFFSUF
    echo "Processing $PREFIX$I"
    echo `python tn_hmm.py -f "$INFILE" -gff "$GFF" > "$INTFILE"`
    echo `python process_genes.py -f "$INTFILE" > "$OUTFILE"`
    echo "$PREFIX$I processed."
    let I++
done
echo "Thank you for using Ben's Amazing TnSeq Analysis Pipline :D"
echo "Enjoy your results, and have a nice day!"
