# 2020 Benjamin J Perry - Attribution-NonCommercial-ShareAlike 4.0 International
# (CC BY-NC-SA 4.0)
# Version: 1.3.0
# Maintainer: Benjamin J Perry
# Email: benjamin.perry@postgrad.otago.ac.nz
# Status: Dev

# Overview:
# R code for handling tnseq wig files and running MANOVA comparison between multiple treatments


readWig() <- function(){
#Reads in a .wig file

}

mergeWig() <- function(){
#Take multipl wig objects and merge into a table

#Option to filter '0'
}

readGFF() <- function(){
#reads in a gff file

}

NormWigTable() <- function(){
#return the clr transformed read counts across insertions sites

}

summarizeFeature() <- function(){
#Take entry in .gff file and returns the normalized read counts and normalize insertion abundance

}

tnseqMANOVA() <- function(){
#computes MANOVA test for feature
#returns unadjusted p-value

}

tnseqFDR() <- function(){
#takes a table of MANOVA p-values and does FDR correction

}
