# __author__ = "Benjamin J Perry"
# __copyright__ = "(CC BY-NC-SA 4.0)"
# __license__ = "MIT License"
# __version__ = "1.0.0"
# __maintainer__ = "Ben Perry"
# __email__ = "ben.perry@agresearch.co.nz"
# __status__ = "alpha"

# Script for generating tradis {sample}.insertionplot.gz files from provided {sample}.bam file

# Usage: ./workflow/scripts/makeInsertionplots.R {input}.bam {input}.bed {wildcards.sample}

# Begin

library(tidyverse) # data manipulation
library(psych) # computing statistics
library(Rsamtools) # extracting sequence length data

make_wig <- function( repliconBed, strand="all", min_count=1, capping=T ){
  
  require(tidyverse)
  require(psych)
  
  if ( length(unique(repliconBed$Chrom)) > 1 ){
    stop("'repliconBed' contains multiple replicons.")
  }
  
  replicon <- as.character(unique(repliconBed$Chrom))
  write( paste("Total Reads in", replicon, ":", length(repliconBed$Chrom), "\n"), stdout())

  #filter the replicon bed file for the indicated strand.
  if ( !(strand %in% c("-", "+", "all")) ){
    stop("'strand' contains invalid character.")
  }
  
  
  write(paste("Filtering ", replicon, " for reads mapped to ", strand, " strand...", sep = ""), stdout())
  if ( strand == "all" ) {
    strandBed <- repliconBed
  } else {
    strandBed <- repliconBed %>% filter( Strand == strand )
  }
  
  
  write(paste(length(strandBed$ChromStart), " reads mapped to ", strand, " strand...", sep = ""), stdout())
  if (strand == "+") {
    # 9 bp offset for tn5 IR element duplication
    wig <- as.data.frame(table(strandBed$ChromStart + 9))
    wig <- wig %>% select("pos" = Var1, "count" = Freq)
    
  } else if (strand == "-") {
    # 9 bp offset for tn5 IR element duplication
    wig <- as.data.frame(table(strandBed$ChromEnd - 9))
    wig <- wig %>% select("pos" = Var1, "count" = Freq)
    
  } else if (strand == "all") {
    # 9 bp offset for tn5 IR element duplication
    wig <- ifelse(strandBed$Strand == "+",
                  strandBed$ChromStart + 9,
                  strandBed$ChromEnd - 9)
    wig <- as.data.frame(table(wig))
    wig <- wig %>% select("pos" = wig, "count" = Freq)
    
  } else {
    stop('incorrect strand value passed.')
  }
  
  write(paste(length(wig$pos), " unique insertions before filtering...", sep = ""), stdout())
  wig$pos <- as.integer(levels(wig$pos))
  
  
  #filter wig by min_count, default is 1
  wig <- wig %>% filter(count >= min_count)
  write(paste(length(wig$pos), " unique insertions after filtering at minimum count threshold of ", min_count, "...", sep = ""), stdout())
  
  
  #Cap the counts at the 99.995% percentile
  if( capping ){
    wig <- cap_wig(wig)
  }
  
  
  write(paste("Unique Insertions with min_count=" , min_count, " and capping=", capping, "...", sep = ""), stdout())
  cat(paste("Summary Statistics for ", replicon, " on ", strand, " strand:\n"), sep = "")
  print(describe(wig[-1]), digits = 2, quote = FALSE)
  
  
  if ( strand == "-" ) {
    wig$count <- as.numeric(wig$count * (-1))
  }
  
  wig <- wig %>% transmute(pos, count = round(count, digits = 0))
  
  
  return(wig)
  
}

cap_wig <- function(wig, SD_threshold=3.5, gMean=NULL, gSDev=NULL){
  
  require(tidyverse)
  wig$logNorm <- log(wig$count)
  
  if ( is.null(gSDev) ) {
    gSDev <- sd(wig$logNorm)
  }
  if ( is.null(gMean) ) {
    gMean <- mean(wig$logNorm)
  }
  
  
  write(paste("Capping counts at ", SD_threshold, " SDs above the geometric mean..."), stdout())
  cap <- gMean + gSDev*SD_threshold
  
  # capping
  wig$logNormCap <- wig$logNorm
  wig <- mutate(wig, logNormCap = ifelse(logNormCap > cap, exp(cap), exp(logNormCap)))
  wig <- wig %>% transmute(pos, "count" = as.numeric(logNormCap))
  
  return(wig)
}

returnBED <- function(BEDFilePath){
  library(tidyverse)
  write(paste("Reading in:", BEDFilePath), stdout())
  bedFile <- read_tsv(BEDFilePath, col_names = F, trim_ws = T, progress = T)
  names(bedFile) <- c("Chrom",
                      "ChromStart",
                      "ChromEnd",
                      "ReadID",
                      "MPQ",
                      "Strand")
  
  return(bedFile)
}

getRepliconLengths <- function( bamFilePath ){
  
  library( Rsamtools )
  # Open .bam file
  BF <- BamFile(file = bamFilePath)
  # Extract replicon names and lengths in dataframe
  replicons <- as.data.frame(seqinfo(BF))
  replicons$replicon <- rownames(replicons)
  #return dataframe
  return( replicons )
  
}

tradis_plot <- function(repliconBed, repliconLength=0, min_count=1, capping=T){
  require(tidyverse)
  require(psych)
  
  
  if ( length( unique( repliconBed$Chrom ) ) > 1 ) {
    stop("'repliconBed' contains multiple replicons.")
  }
  
  
  replicon <- as.character(unique(repliconBed$Chrom))
  write(paste("Total Reads in", replicon, ":", length(repliconBed$Chrom)), "\n", stdout())

  # make plus strand with tn5 offset
  pBed <- repliconBed %>% filter(Strand == "+")
  pWig <- as.data.frame(table(pBed$ChromStart + 10))
  rm(pBed)
  
  
  pWig <- pWig %>% filter(Freq >= min_count) %>% 
    mutate(
      "pos" = as.integer(as.character(Var1)),
      "count" = Freq,
      Var1 = NULL,
      Freq = NULL
    )
  
  # make minus strand with tn5 offset
  mBed <- repliconBed %>% filter(Strand == "-")
  mWig <- as.data.frame(table(mBed$ChromEnd - 10))
  rm(mBed)
  
  mWig <- mWig %>% filter(Freq >= min_count) %>% 
    mutate(
      "pos" = as.integer(as.character(Var1)),
      "count" = Freq,
      Var1 = NULL,
      Freq = NULL
    )
  
  # make all wig by merging pWig and mWig for capping parameters
  allWig <- bind_rows(mWig, pWig)
  describe(allWig$count)
  
  # capping parameters for stranded wigs
  gSDev <- sd(log(allWig$count))
  gMean <- mean(log(allWig$count))
  
  # capping the stranded wigs
  pcWig <- cap_wig(wig = pWig, gMean = gMean, gSDev = gSDev)
  mcWig <- cap_wig(wig = mWig, gMean = gMean, gSDev = gSDev)
  
  
  ### Merge and return the TraDIS insetion tack
  tradisPlot <- data.frame( pos = as.integer(seq(from = 1, to = repliconLength, by = 1)) )
  
  tradisPlot <- left_join(x = tradisPlot, y = pcWig, "pos")
  tradisPlot <- left_join(x = tradisPlot, y = mcWig, "pos")
  
  
  # Clean up column names and NA values = 0
  tradisPlot <- tradisPlot %>% mutate(
    "plus" = as.integer(replace_na(count.x, replace = 0)),
    "minus" = as.integer(replace_na(count.y, replace = 0)),
    count.x = NULL,
    count.y = NULL
  )
  
  
  return(tradisPlot)
}

write_insertionplots <- function(sampleID, repliconID, repliconLength, bedFile, outfile){
  
  library(tidyselect)
  
  repliconBed <- bedFile %>% dplyr::filter( Chrom == repliconID )
  
  tradisPlot <- tradis_plot( repliconBed = repliconBed, min_count = 1, repliconLength = repliconLength, capping = T )
  
  write_delim(tradisPlot[-1], outfile, delim = " ", col_names = F)
  
  return(write(paste("Tradis .insert_site_plot.gz written to", outfile, "\n", sep = " "), stdout()))
}

# ### Parse Arguments
args <- commandArgs( trailingOnly = TRUE )

if ( length(args) == 0 ){
    stop("Error: No command line arguments given.")
}

BAMFilePath <- args[1]
if (  !grepl(".bam", BAMFilePath) ){
    stop("Error: bam file path invalid; No '.bam'")
}
if (  !file.exists(BAMFilePath) ){
    stop("Error: .bam file does not exsist")
}

BEDFilePath <- args[2]
if ( !grepl(".bed", BEDFilePath) ){
    stop("Error: bam file path invalid; No '.bam'")
}
if (  !file.exists(BEDFilePath) ){
    stop("Error: .bed file does not exsist")
}

SAMPLE <- args[3]
if ( length(SAMPLE) == 0 ){
    stop("Error: No sample given.")
}

# ### Read in .bed file for processing

bedfile <- returnBED(BEDFilePath)


# ## Parse the bam for replicons and seqlengths

replicons <- getRepliconLengths( BAMFilePath )

replicons$sample <- SAMPLE


# ### make TraDIS insertion plot

# Iterate over replicons and make insert_site_plot.gz files
outdirInsertions <- "output/05_tradis_plots/"
if ( !dir.exists(outdirInsertions) ){
  dir.create(path = outdirInsertions, recursive = TRUE)
} else {
  write(paste(outdirInsertions, "output directory already exsists. Proceeding.", sep = " "), stdout())
}

for(i in 1:nrow(replicons)) { 
  
  # Prepare the output file handle
  tPlotOut <- paste( outdirInsertions, replicons[i, "sample"], ".", replicons[i, "replicon"], ".insert_site_plot.gz", sep = "" ) 
  
  file.create( tPlotOut )
  
  # make the replicon insert_site_plot
  write_insertionplots(sampleID = replicons[i, "sample"], 
                       repliconID = replicons[i, "replicon"], 
                       repliconLength = replicons[i, "seqlengths"], 
                       bedFile = bedfile,
                       outfile = tPlotOut)
}


# ### Making wig files

outdirWigs <- "output/06_wig_files/"
if ( !dir.exists(outdirWigs) ){
  dir.create(path = outdirWigs, recursive = TRUE)
} else {
  write(paste(outdirWigs, "output directory already exsists. Proceeding.", sep = " "), stdout())
}

pWigOut <- paste(outdirWigs, replicons[i, "sample"], ".plus.wig", sep = "")
mWigOut <- paste(outdirWigs, replicons[i, "sample"], ".minus.wig", sep = "")
filtWigOut <- paste(outdirWigs, replicons[i, "sample"], ".filt.1.wig", sep = "")

file.create(pWigOut)
file.create(mWigOut)
file.create(filtWigOut)


for(i in 1:nrow(replicons)) { 
  
  replicon <- replicons[i, "replicon"]
  
  headerLine <- paste("variableStep chrom=", replicon, sep = '')
  
  repliconBed <- bedfile %>% dplyr::filter(Chrom == replicon)
  
  # plus strand wig file
  plusWig <- make_wig(repliconBed = repliconBed, 
                      strand = "+", 
                      min_count = 1, 
                      capping = T)
  

  write(paste("Printing .wig file:", pWigOut, "\n"), stdout())
  
  cat(headerLine, file = pWigOut, sep = "\n", append = T)
  
  write_delim(x = plusWig,
              file = pWigOut,
              delim = " ",
              col_names = F,
              append = T)

  # minus strand wig file
  minusWig <-make_wig(repliconBed = repliconBed,
                      strand = "-",
                      min_count = 1,
                      capping = T)
  write(paste("Printing .wig file:", mWigOut, "\n"), stdout())
  
  cat(headerLine,
      file = mWigOut,
      sep = "\n",
      append = T)
  
  write_delim(x = minusWig,
              file = mWigOut,
              delim = " ",
              col_names = F,
              append = T)

  # all wig file
  allWig <- make_wig(repliconBed = repliconBed,
                     strand = "all",
                     min_count = 1,
                     capping = T)
  write(paste("Printing .wig file:", filtWigOut, "\n"), stdout())
  
  cat(headerLine,
      file = filtWigOut,
      sep = "\n",
      append = T)
  
  write_delim(x = allWig,
              file = filtWigOut,
              delim = " ",
              col_names = F,
              append = T)
  write("\n\n\n", stdout())
  
}
