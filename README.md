# Scripts and Archive for Tnseq Analysis

## Description

This repo contains scripts necessary to analyze and visualize TnSeq data.

## scripts/RefGenTA.py

Utility script for computing 'TA' motif reference wig files visualization.
Also a helper script for WigScripts.py when generating .wif tracks for treatment .bed files.

## scripts/WigScripts.py

Utility script for computing treatment .wig tracks. Uses reference .fasta file used for alignment of tn-tags,
a .bed file of the aligned tn-tags, and a for an output .wig track to be printed to.

## Usage

Both RegGenTA.py and WigScripts.py have command line usage. use '-h' flag for help menus.