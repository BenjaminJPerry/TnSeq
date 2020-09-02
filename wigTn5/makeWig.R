library(tidyverse)
library(psych)


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
  
  strandBed <- repliconBed %>% filter(Strand == strand)
  
  if (strand == "+") {
    wig <- as.data.frame(table(strandBed$ChromStart + 10))
  } else if (strand == "-") {
    wig <- as.data.frame(table(repliconBed$ChromEnd - 10))
  } else if (strand == "all") {
    obsAll <- ifelse(repliconBed$Strand == "+",
                     repliconBed$ChromStart + 10,
                     repliconBed$ChromEnd - 10)
    wig <- as.data.frame(table(obsAll))
  } else {
    stop('incorrect strand value passed.')
  }
  wig$Var1 <- as.integer(levels(wig$Var1))
  wig <- wig %>% select("pos" = Var1, "count" = Freq)
  
  #filter wig by min_count, default is 3
  wig <- wig %>% filter(count >= min_count)
  
  #Cap the counts at the 99.995% percentile
  if(capping){
    wig <- cap_wig(wig)
  }
  
  print(paste("Unique Insertions with min_count=" , min_count, " and capping=", capping, ": ", length(wig$count), sep = ""))
  cat("Summary Statistics:\n")
  print(describe(wig[-1]), digits = 2)
  
  if (strand == "-") {
  wig$logNormCap <- wig$logNormCap * (-1)
  }
  
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
  
  cap <- gMean + gSDev*SD_threshold
  
  # capping
  wig$logNormCap <- wig$logNorm
  wig <- mutate(wig, logNormCap = ifelse(logNormCap > cap, exp(cap), exp(logNormCap)))
  wig <- wig %>% select(pos, "count" = logNormCap)
  
  return(wig)
}

#TODO: make TraDIS insertion plot
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
  pWig <- pWig %>% filter(Freq >= min_count) %>% mutate("pos" = as.integer(as.character(Var1)), "count" = Freq, Var1 = NULL, Freq = NULL)
  
  # make minus strand
  mBed <- repliconBed %>% filter(Strand == "-")
  mWig <- as.data.frame(table(mBed$ChromEnd - 10))
  rm(mBed)
  mWig <- mWig %>% filter(Freq >= min_count) %>% mutate("pos" = as.integer(as.character(Var1)), "count" = Freq, Var1 = NULL, Freq = NULL)
  
  # make all wig
  allWig <- bind_rows(mWig, pWig)
  describe(allWig$count)
  gSDev <- sd(log(allWig$count))
  gMean <- mean(log(allWig$count))
  
  pcWig <- cap_wig(wig = pWig, gMean = gMean, gSDev = gSDev)
  mcWig <- cap_wig(wig = mWig, gMean = gMean, gSDev = gSDev)
  
  ### Merge and return the TraDIS insetion tack
}



for (seqLib in list.dirs(recursive = F)) {
        #Prepareing the directory tree and files paths for the loop.
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
        
        
        #TODO: Revise this sectoin of the loop once the functions have been written.
        pWigOut <- paste(wigOutPath, "/", seqRoot, ".plus.wig", sep = "")
        mWigOut <- paste(wigOutPath, "/", seqRoot, ".minus.wig", sep = "")
        allWigOut <- paste(wigOutPath, "/", seqRoot, ".all.wig", sep = "")
        filtWigOut <- paste(wigOutPath, "/", seqRoot, ".filt.1.wig", sep = "")
        file.create(pWigOut)
        file.create(mWigOut)
        file.create(allWigOut)
        file.create(filtWigOut)
        
        #Nested loop for unique replicon in the bedFile.  
        for (entry in unique(bedFile$Chrom)) {
                replicon <- as.character(entry)
                headerLine <- paste("variableStep chrom=", replicon, sep = '')
                
                repliconBed <- bedFile %>% filter(Chrom == replicon)
                

                print(paste("Total Reads in", replicon, ":", length(repliconBed$Chrom)))
                cat("\n")
                PlusStrand <- repliconBed %>% filter(Strand == "+")
                MinusStrand <- repliconBed %>% filter(Strand == "-")

                #Plus Strand Wig Stats and Wig File
                PlusWig <- as.data.frame(table(PlusStrand$ChromStart + 10))
                print(paste("Printing .wig file:", pWigOut))
                print(paste("Unique Insertions Plus Strand: ", length(PlusWig$Freq)))
                cat("Summary Statistics:\n")
                print(describe(PlusWig$Freq, trim = 0.05, IQR = T))
                PlusWig$Var1 <- as.integer(levels(PlusWig$Var1))

                cat(headerLine,
                    file = pWigOut,
                    sep = "\n",
                    append = T
                    )

                write_delim(
                        x = PlusWig,
                        path = pWigOut,
                        delim = " ",
                        col_names = F,
                        append = T
                )
                cat("\n\n")


                #Minus Strand Wig Stats and Wig File
                MinusWig <- as.data.frame(table(MinusStrand$ChromEnd - 10))
                MinusWigforDescribe <- MinusWig
                MinusWig$Freq <- MinusWig$Freq * (-1)
                print(paste("Printing .wig file:", mWigOut))
                print(paste(
                        "Unique Insertions Minus Strand: ",
                        length(MinusWigforDescribe$Freq)
                ))
                cat("Summary Statistics:\n")
                print(describe(
                        MinusWigforDescribe$Freq,
                        trim = 0.05,
                        IQR = T
                ))
                MinusWig$Var1 <- as.integer(levels(MinusWig$Var1))


                cat(headerLine,
                    file = mWigOut,
                    sep = "\n",
                    append = T
                )


                write_delim(
                        x = MinusWig,
                        path = mWigOut,
                        delim = " ",
                        col_names = F,
                        append = T
                )
                cat("\n\n")


                ObsAll <-ifelse(repliconBed$Strand == "+",
                                repliconBed$ChromStart + 10,
                                repliconBed$ChromEnd - 10) # Combined offsets
                AllWig <- as.data.frame(table(ObsAll))

                print(paste("Printing .wig file:", allWigOut))
                print(paste("Unique Insertions In Total: ", length(AllWig$Freq)))
                cat("Summary Statistics:\n")
                print(describe(AllWig$Freq, trim = 0.05, IQR = T))
                AllWig$ObsAll <-
                        as.integer(levels(AllWig$ObsAll))

                cat(headerLine,
                    file = allWigOut,
                    sep = "\n",
                    append = T
                    )

                write_delim(
                        x = AllWig,
                        path = allWigOut,
                        delim = " ",
                        col_names = F,
                        append = T
                )
                cat("\n\n")


                FiltAll <- AllWig %>% filter(Freq > 1) #remove singletons
                print(paste("Printing .wig file:", filtWigOut))
                print(paste("Unique Insertions In Total: ", length(FiltAll$Freq)))
                cat("Summary Statistics:\n")
                print(describe(FiltAll$Freq, trim = 0.05, IQR = T))

                cat(headerLine,
                    file = filtWigOut,
                    sep = "\n",
                    append = T
                    )

                write_delim(
                        x = FiltAll,
                        path = filtWigOut,
                        delim = " ",
                        col_names = F,
                        append = T
                )
                cat("\n\n")

                #Clean large variable from memory
                # rm(bedFile)
                # rm(PlusStrand)
                # rm(MinusStrand)
                # rm(PlusWig)
                # rm(MinusWig)
                # rm(ObsAll)
                # rm(AllWig)
                # rm(FiltAll)
                # rm(FiltAll2)
        }
}
