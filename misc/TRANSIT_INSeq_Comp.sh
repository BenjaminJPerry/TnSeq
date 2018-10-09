#!/bin/bash

###Transit Analysis Automation Pipeline
PROTABLE="/home/benji/Desktop/Transit_Analysis/RLV3841.prot_table"
INPUTBASE="/home/benji/Desktop/Transit_Analysis/0hrSpermosphere/AM236080"
OXFINSEQ="/home/benji/Desktop/Transit_Analysis"

INPUT1="$INPUTBASE.0hrSPOS.1.wig"
INPUT2="$INPUTBASE.0hrSPOS.2.wig"
INPUT3="$INPUTBASE.0hrSPOS.3.wig"
INPUT="$INPUT1,$INPUT2,$INPUT3" ### 0hrSPOS input wigs

REPLICON="/AM236080"

OUTDIR="/0hrSPOS_Resampling_Analysis"

### Resampling Analysis TYVMM vs 0hrSPOS
TREAT=".vmmty"
SUBDIR="/VMMTYInput"
CONTROL1="$OXFINSEQ$SUBDIR$REPLICON$TREAT.1.wig"
CONTROL2="$OXFINSEQ$REPLICON$TREAT.2.wig"
CONTROL3="$OXFINSEQ$REPLICON$TREAT.3.wig"
CONTROL="$CONTROL1,$CONTROL1,$CONTROL1"
OUTFILE="$OXFINSEQ$OUTDIR$REPLICON$TREAT.0hrSPOS.dat"

python /home/benji/Research_Programs/transit/src/transit.py resampling $PROTABLE $CONTROL $INPUT $OUTFILE -s 10000 -H -N TTR -L -iN 0 -iC 0

### Resampling Analysis WaterAgar vs 0hrSPOS
TREAT=".wateragar"
SUBDIR="/48hrWaterAgar"
CONTROL1="$OXFINSEQ$SUBDIR$REPLICON$TREAT.1.wig"
CONTROL2="$OXFINSEQ$REPLICON$TREAT.2.wig"
CONTROL3="$OXFINSEQ$REPLICON$TREAT.3.wig"
CONTROL="$CONTROL1,$CONTROL1,$CONTROL1"
OUTFILE="$OXFINSEQ$OUTDIR$REPLICON$TREAT.0hrSPOS.dat"

python /home/benji/Research_Programs/transit/src/transit.py resampling $PROTABLE $CONTROL $INPUT $OUTFILE -s 10000 -H -N TTR -L -iN 0 -iC 0

### Resampling Analysis IVR vs 0hrSPOS
TREAT=".IVR"
SUBDIR="/InVitroRadicle"
CONTROL1="$OXFINSEQ$SUBDIR$REPLICON$TREAT.1.wig"
CONTROL2="$OXFINSEQ$REPLICON$TREAT.2.wig"
CONTROL3="$OXFINSEQ$REPLICON$TREAT.3.wig"
CONTROL="$CONTROL1,$CONTROL1,$CONTROL1"
OUTFILE="$OXFINSEQ$OUTDIR$REPLICON$TREAT.0hrSPOS.dat"

python /home/benji/Research_Programs/transit/src/transit.py resampling $PROTABLE $CONTROL $INPUT $OUTFILE -s 10000 -H -N TTR -L -iN 0 -iC 0


### Resampling Analysis 48hrSPOS vs 0hrSPOS
TREAT=".48hrSPOS"
SUBDIR="/48hrSpermosphere"
CONTROL1="$OXFINSEQ$SUBDIR$REPLICON$TREAT.1.wig"
CONTROL2="$OXFINSEQ$REPLICON$TREAT.2.wig"
CONTROL3="$OXFINSEQ$REPLICON$TREAT.3.wig"
CONTROL="$CONTROL1,$CONTROL1,$CONTROL1"
OUTFILE="$OXFINSEQ$OUTDIR$REPLICON$TREAT.0hrSPOS.dat"

python /home/benji/Research_Programs/transit/src/transit.py resampling $PROTABLE $CONTROL $INPUT $OUTFILE -s 10000 -H -N TTR -L -iN 0 -iC 0
