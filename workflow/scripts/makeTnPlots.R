# 2020 Benjamin J Perry - Attribution-NonCommercial-ShareAlike 4.0 International
# (CC BY-NC-SA 4.0)
# Version: 2.0.0
# Maintainer: Benjamin J Perry
# Email: benjamin.perry@postgrad.otago.ac.nz
# Status: Functional

library(tidyverse)
library(psych)
options(width=120)

# wig generating function
make_wig <- function(repliconBed, strand="all", min_count=3, capping=T){
  require(tidyverse)
  require(psych)
  if (length(unique(repliconBed$Chrom)) > 1) stop("'repliconBed' contains multiple replicons.")
  replicon <- as.character(unique(repliconBed$Chrom))

  print(paste("Total Reads in", replicon, ":", length(repliconBed$Chrom)))
  cat("\n")

  #filter the replicon bedList for the indicated strand.
  if (!(strand %in% c("-", "+", "all"))) stop("'strand' contains invalid character.")

  print(paste("Filtering ", replicon, " for reads mapped to ", strand, "...", sep = ""))

  if (strand == "all") {
    strandBed <- repliconBed
  } else {
    strandBed <- repliconBed %>% filter(Strand == strand)
  }

  print(paste(length(strandBed$ChromStart), " reads mapped to ", strand, "...", sep = ""))

  if (strand == "+") {
    wig <- as.data.frame(table(strandBed$ChromStart + 9))
    wig <- wig %>% select("pos" = Var1, "count" = Freq)
  } else if (strand == "-") {
    wig <- as.data.frame(table(strandBed$ChromEnd - 9))
    wig <- wig %>% select("pos" = Var1, "count" = Freq)
  } else if (strand == "all") {
    wig <- ifelse(strandBed$Strand == "+",
                  strandBed$ChromStart + 9,
                  strandBed$ChromEnd - 9)
    wig <- as.data.frame(table(wig))
    wig <- wig %>% select("pos" = wig, "count" = Freq)
  } else {
    stop('incorrect strand value passed.')
  }
  print(paste(length(wig$pos), " unique insertions before filtering...", sep = ""))
  wig$pos <- as.integer(levels(wig$pos))

  #filter wig by min_count, default is 3
  wig <- wig %>% filter(count >= min_count)
  print(paste(length(wig$pos), " unique insertions after filtering at minimum count threshold of ", min_count, "...", sep = ""))

  #Cap the counts at the 99.995% percentile
  if(capping){
    wig <- cap_wig(wig)
  }

  print(paste("Unique Insertions with min_count=" , min_count, " and capping=", capping, ": ", length(wig$count), sep = ""))
  cat(paste("Summary Statistics for ", replicon, " on ", strand, " strand:\n"), sep = "")
  print(describe(wig[-1]), digits = 2)

  if (strand == "-") {
  wig$count <- as.numeric(wig$count * (-1))
  }

  wig <- wig %>% transmute(pos, count = round(count, digits = 0))
  return(wig)

}
# capping function for insertion track
cap_wig <- function(wig, SD_threshold=3.5, gMean=NULL, gSDev=NULL){
  require(tidyverse)
  wig$logNorm <- log(wig$count)

  if (is.null(gSDev)) {
    gSDev <- sd(wig$logNorm)
  }

  if (is.null(gMean)) {
    gMean <- mean(wig$logNorm)
  }
  print(paste("Capping counts at ", SD_threshold, " SDs above the mean..."))
  cap <- gMean + gSDev*SD_threshold

  # capping
  wig$logNormCap <- wig$logNorm
  wig <- mutate(wig, logNormCap = ifelse(logNormCap > cap, exp(cap), exp(logNormCap)))
  wig <- wig %>% transmute(pos, "count" = as.numeric(logNormCap))

  return(wig)
}

# make TraDIS insertion plot
tradis_plot <- function(repliconBed, repliconLength=6530403, min_count=3, capping=T){
  # replicon length is currently hardcoded for R7A
  require(tidyverse)
  require(psych)
  if (length(unique(repliconBed$Chrom)) > 1) stop("'repliconBed' contains multiple replicons.")
  replicon <- as.character(unique(repliconBed$Chrom))

  print(paste("Total Reads in", replicon, ":", length(repliconBed$Chrom)))
  cat("\n")

  # make plus strand
  pBed <- repliconBed %>% filter(Strand == "+")
  pWig <- as.data.frame(table(pBed$ChromStart + 10))
  rm(pBed)
  pWig <-
    pWig %>% filter(Freq >= min_count) %>% mutate(
      "pos" = as.integer(as.character(Var1)),
      "count" = Freq,
      Var1 = NULL,
      Freq = NULL
    )

  # make minus strand
  mBed <- repliconBed %>% filter(Strand == "-")
  mWig <- as.data.frame(table(mBed$ChromEnd - 10))
  rm(mBed)
  mWig <-
    mWig %>% filter(Freq >= min_count) %>% mutate(
      "pos" = as.integer(as.character(Var1)),
      "count" = Freq,
      Var1 = NULL,
      Freq = NULL
    )

  # make all wig
  allWig <- bind_rows(mWig, pWig)
  describe(allWig$count)
  gSDev <- sd(log(allWig$count))
  gMean <- mean(log(allWig$count))

  pcWig <- cap_wig(wig = pWig, gMean = gMean, gSDev = gSDev)
  mcWig <- cap_wig(wig = mWig, gMean = gMean, gSDev = gSDev)

  ### Merge and return the TraDIS insetion tack
  tradisPlot <- data.frame(
    pos = as.integer(seq(from = 1, to = repliconLength, by = 1))
  )
  tradisPlot <- left_join(x = tradisPlot, y = pcWig, "pos")
  tradisPlot <- left_join(x = tradisPlot, y = mcWig, "pos")
  # Clean up column names and NA values = 0
  tradisPlot <-
    tradisPlot %>% mutate(
      "plus" = as.integer(replace_na(count.x, replace = 0)),
      "minus" = as.integer(replace_na(count.y, replace = 0)),
      count.x = NULL,
      count.y = NULL
    )

  return(tradisPlot)
}

#Loop for each directory and nested loop for each replicon in bedFile
for (seqLib in list.dirs(recursive = F)) {
        #Preparing the directory tree and files paths for the loop.
        cat("\n\n\n")
        seqRoot <- str_split(seqLib, "/", simplify = T)[2]
        print(paste("Processing Sample: ", seqRoot))
        bedIN <- paste(seqLib, "/alignment/", seqRoot, ".bed", sep = '')
        print(paste("Reading in:", bedIN), quote = F)
        bedFile <- read_tsv(bedIN, col_names = F, trim_ws = T, progress = T)
        names(bedFile) <- c("Chrom",
                            "ChromStart",
                            "ChromEnd",
                            "ReadID",
                            "MPQ",
                            "Strand")

        wigOutPath <- paste(seqLib, "/wig", sep = '')
        dir.create(path = wigOutPath, showWarnings = T)

        pWigOut <- paste(wigOutPath, "/", seqRoot, ".plus.wig", sep = "")
        mWigOut <- paste(wigOutPath, "/", seqRoot, ".minus.wig", sep = "")
        filtWigOut <- paste(wigOutPath, "/", seqRoot, ".filt.3.wig", sep = "")
        file.create(pWigOut)
        file.create(mWigOut)
        file.create(filtWigOut)

        tPlotOutPath <- paste(seqLib, "/tradis", sep = '')
        dir.create(path = tPlotOutPath, showWarnings = T)
        tPlotOut <- paste(tPlotOutPath, "/", seqRoot, ".tradis_insertion_plot.gz", sep = "")
        file.create(tPlotOut)

        #Nested loop for unique replicon in the bedFile.
        for (entry in unique(bedFile$Chrom)) {
                replicon <- as.character(entry)
                headerLine <- paste("variableStep chrom=", replicon, sep = '')
                repliconBed <- bedFile %>% filter(Chrom == replicon)

                # plus strand wig file
                plusWig <- make_wig(repliconBed = repliconBed, strand = "+", min_count = 3, capping = T)
                cat("\n")
                print(paste("Printing .wig file:", pWigOut))
                cat(headerLine,
                    file = pWigOut,
                    sep = "\n",
                    append = T
                    )

                write_delim(
                        x = plusWig,
                        path = pWigOut,
                        delim = " ",
                        col_names = F,
                        append = T
                )
                cat("\n")

                # minus strand wig file
                minusWig <- make_wig(repliconBed = repliconBed, strand = "-", min_count = 3, capping = T)
                cat("\n")
                print(paste("Printing .wig file:", mWigOut))
                cat(headerLine,
                    file = mWigOut,
                    sep = "\n",
                    append = T
                )

                write_delim(
                  x = minusWig,
                  path = mWigOut,
                  delim = " ",
                  col_names = F,
                  append = T
                )
                cat("\n")

                # all wig file
                allWig <- make_wig(repliconBed = repliconBed, strand = "all", min_count = 3, capping = T)
                print(paste("Printing .wig file:", filtWigOut))

                cat(headerLine,
                    file = filtWigOut,
                    sep = "\n",
                    append = T
                    )

                write_delim(
                        x = allWig,
                        path = filtWigOut,
                        delim = " ",
                        col_names = F,
                        append = T
                )
                cat("\n")

                # tradis insertion plot
                tradisPlot <- tradis_plot(repliconBed = repliconBed, min_count = 3, capping = T)
                write_delim(tradisPlot[-1], tPlotOut, delim = " ", col_names = F)

        }
}
