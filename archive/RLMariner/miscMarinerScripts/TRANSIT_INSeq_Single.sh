#!/bin/bash

###Transit Analysis Automation Pipeline
PROTABLE="/home/benji/Desktop/Transit_Analysis/RLV3841.prot_table"
printf "Define Sample Name Base: "
read NAMEBASE

WORKINGDIR=$(pwd)

INPUT1="$NAMEBASE.1.wig"
INPUT2="$NAMEBASE.2.wig"
INPUT3="$NAMEBASE.3.wig"
INPUT="$INPUT1,$INPUT2,$INPUT3"


### HMM Analysis with SUM
OUTPUT="$(pwd)/$NAMEBASE.HMM.SUM.dat"
python /home/benji/Research_Programs/transit/src/transit.py hmm $PROTABLE $INPUT $OUTPUT -r Sum -iN 0 -iC 0

### HMM Analysis with MEAN
OUTPUT="$(pwd)/$NAMEBASE.HMM.MEAN.dat"
python /home/benji/Research_Programs/transit/src/transit.py hmm $PROTABLE $INPUT $OUTPUT -r Mean -iN 0 -iC 0

### gumbel
#gumbel using Sum, Min Read from 1-5
OUTPUT="$(pwd)/$NAMEBASE.gumbel.m1.SUM.dat"
python /home/benji/Research_Programs/transit/src/transit.py gumbel $PROTABLE $INPUT $OUTPUT -s 10000 -m 1 -b 500 -t 1 -r Sum -iN 0 -iC 0
OUTPUT="$(pwd)/$NAMEBASE.gumbel.m2.SUM.dat"
python /home/benji/Research_Programs/transit/src/transit.py gumbel $PROTABLE $INPUT $OUTPUT -s 10000 -m 2 -b 500 -t 1 -r Sum -iN 0 -iC 0
OUTPUT="$(pwd)/$NAMEBASE.gumbel.m3.SUM.dat"
python /home/benji/Research_Programs/transit/src/transit.py gumbel $PROTABLE $INPUT $OUTPUT -s 10000 -m 3 -b 500 -t 1 -r Sum -iN 0 -iC 0
OUTPUT="$(pwd)/$NAMEBASE.gumbel.m4.SUM.dat"
python /home/benji/Research_Programs/transit/src/transit.py gumbel $PROTABLE $INPUT $OUTPUT -s 10000 -m 4 -b 500 -t 1 -r Sum -iN 0 -iC 0
OUTPUT="$(pwd)/$NAMEBASE.gumbel.m5.SUM.dat"
python /home/benji/Research_Programs/transit/src/transit.py gumbel $PROTABLE $INPUT $OUTPUT -s 10000 -m 5 -b 500 -t 1 -r Sum -iN 0 -iC 0

#gumbel using Mean, Min Read 1-5
OUTPUT="$(pwd)/$NAMEBASE.gumbel.m1.MEAN.dat"
python /home/benji/Research_Programs/transit/src/transit.py gumbel $PROTABLE $INPUT $OUTPUT -s 10000 -m 1 -b 500 -t 1 -r Mean -iN 0 -iC 0
OUTPUT="$(pwd)/$NAMEBASE.gumbel.m2.MEAN.dat"
python /home/benji/Research_Programs/transit/src/transit.py gumbel $PROTABLE $INPUT $OUTPUT -s 10000 -m 2 -b 500 -t 1 -r Mean -iN 0 -iC 0
OUTPUT="$(pwd)/$NAMEBASE.gumbel.m3.MEAN.dat"
python /home/benji/Research_Programs/transit/src/transit.py gumbel $PROTABLE $INPUT $OUTPUT -s 10000 -m 3 -b 500 -t 1 -r Mean -iN 0 -iC 0
OUTPUT="$(pwd)/$NAMEBASE.gumbel.m4.MEAN.dat"
python /home/benji/Research_Programs/transit/src/transit.py gumbel $PROTABLE $INPUT $OUTPUT -s 10000 -m 4 -b 500 -t 1 -r Mean -iN 0 -iC 0
OUTPUT="$(pwd)/$NAMEBASE.gumbel.m5.MEAN.dat"
python /home/benji/Research_Programs/transit/src/transit.py gumbel $PROTABLE $INPUT $OUTPUT -s 10000 -m 5 -b 500 -t 1 -r Mean -iN 0 -iC 0

exit
