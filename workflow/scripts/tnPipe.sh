# 2020 Benjamin J Perry - Attribution-NonCommercial-ShareAlike 4.0 International
# (CC BY-NC-SA 4.0)
# Version: 1.0.0
# Maintainer: Benjamin J Perry
# Email: benjamin.perry@postgrad.otago.ac.nz
# Status: Development

#!/bin/sh

~/projects/TnSeq/wigTn5/makeBed.sh

Rscript --verbose ~/projects/TnSeq/wigTn5/makeTnPlots.R

### TraDIS Analysis of insertion plots
R7AEMBL=~/ref/R7A_20-4-202000000000.current.embl
for i in $(ls);
do
  python ~/projects/TnSeq/tnScripts/wigScripts.py -F ~/ref/R7A_20-4-202000000000.current.fasta -B "$i"/alignment/"$i".bed -O "$i"/wig/"$i".tn5.wig -Tn5

  cd $i/tradis
  TRADISPLOT=$i.tradis_insertion_plot.gz

  tradis_gene_insert_sites $R7AEMBL -trim5 0.05 -trim3 0.1 $TRADISPLOT

  TRADISINS=$i.tradis_insertion_plot.tradis_gene_insert_sites.csv
  #TRADISICE=$i.tradis_insertion_plot.tradis_gene_insert_sites.ICE.csv
  TRADISICEFILT=$i.tradis_insertion_plot.tradis_gene_insert_sites.ICE.Filt.csv
  ICELOCI=~/ref/ICESym.locus.tags.txt


  #grep -f $ICELOCI $i.tradis_insertion_plot.tradis_gene_insert_sites.csv > $TRADISICE
  grep -vf $ICELOCI $i.tradis_insertion_plot.tradis_gene_insert_sites.csv > $TRADISICEFILT

  tradis_essentiality.R $TRADISINS
  tradis_essentiality.R $TRADISICEFILT

  cd ../../

done
exit
