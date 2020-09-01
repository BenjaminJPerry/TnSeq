library(tidyverse)
library(psych)

for (seqLib in list.dirs(recursive = F)) {
        cat("\n\n\n")
        seqRoot <- str_split(seqLib, "/", simplify = T)[2]
        print(paste("Processing Sample: ", seqRoot))
        bedIN <- paste(seqLib, "/alignment/", seqRoot, ".bed", sep = '')
        print(paste("Reading in:", bedIN), quote = F)
        BedFile <- read_tsv(bedIN, col_names = F, trim_ws = T, progress = T)
        names(BedFile) <- c("Chrom",
                            "ChromStart",
                            "ChromEnd",
                            "ReadID",
                            "MPQ",
                            "Strand")
        
        wigOutPath <- paste(seqLib, "/wig", sep = '')
        dir.create(path = wigOutPath, showWarnings = T)
        
        PWigOut <- paste(wigOutPath, "/", seqRoot, ".plus.wig", sep = "")
        MWigOut <- paste(wigOutPath, "/", seqRoot, ".minus.wig", sep = "")
        AllWigOut <- paste(wigOutPath, "/", seqRoot, ".all.wig", sep = "")
        FiltWigOut <- paste(wigOutPath, "/", seqRoot, ".filt.1.wig", sep = "")
        file.create(PWigOut)
        file.create(MWigOut)
        file.create(AllWigOut)
        file.create(FiltWigOut)
        
        for (entry in unique(BedFile$Chrom)) {
                replicon <- as.character(entry)
                repliconBed <- BedFile %>% filter(Chrom == replicon)
                headerLine <- paste("variableStep chrom=", replicon, sep = '')
                
                print(paste("Total Reads in", replicon, ":", length(repliconBed$Chrom)))
                cat("\n")
                PlusStrand <- repliconBed %>% filter(Strand == "+")
                MinusStrand <- repliconBed %>% filter(Strand == "-")
                
                #Plus Strand Wig Stats and Wig File
                PlusWig <- as.data.frame(table(PlusStrand$ChromStart + 10))
                print(paste("Printing .wig file:", PWigOut))
                print(paste("Unique Insertions Plus Strand: ", length(PlusWig$Freq)))
                cat("Summary Statistics:\n")
                print(describe(PlusWig$Freq, trim = 0.05, IQR = T))
                PlusWig$Var1 <- as.integer(levels(PlusWig$Var1))
                
                cat(headerLine,
                    file = PWigOut,
                    sep = "\n",
                    append = T
                    )
                
                write_delim(
                        x = PlusWig,
                        path = PWigOut,
                        delim = " ",
                        col_names = F,
                        append = T
                )
                cat("\n\n")
                
                
                #Minus Strand Wig Stats and Wig File
                MinusWig <- as.data.frame(table(MinusStrand$ChromEnd - 10))
                MinusWigforDescribe <- MinusWig
                MinusWig$Freq <- MinusWig$Freq * (-1)
                print(paste("Printing .wig file:", MWigOut))
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
                    file = MWigOut,
                    sep = "\n",
                    append = T
                )
                
                
                write_delim(
                        x = MinusWig,
                        path = MWigOut,
                        delim = " ",
                        col_names = F,
                        append = T
                )
                cat("\n\n")
                
                
                ObsAll <-ifelse(repliconBed$Strand == "+",
                                repliconBed$ChromStart + 10,
                                repliconBed$ChromEnd - 10) # Combined offsets
                AllWig <- as.data.frame(table(ObsAll))
                
                print(paste("Printing .wig file:", AllWigOut))
                print(paste("Unique Insertions In Total: ", length(AllWig$Freq)))
                cat("Summary Statistics:\n")
                print(describe(AllWig$Freq, trim = 0.05, IQR = T))
                AllWig$ObsAll <-
                        as.integer(levels(AllWig$ObsAll))
                
                cat(headerLine,
                    file = AllWigOut,
                    sep = "\n",
                    append = T
                    )
                
                write_delim(
                        x = AllWig,
                        path = AllWigOut,
                        delim = " ",
                        col_names = F,
                        append = T
                )
                cat("\n\n")
                
                
                FiltAll <- AllWig %>% filter(Freq > 1) #remove singletons
                print(paste("Printing .wig file:", FiltWigOut))
                print(paste("Unique Insertions In Total: ", length(FiltAll$Freq)))
                cat("Summary Statistics:\n")
                print(describe(FiltAll$Freq, trim = 0.05, IQR = T))
                
                cat(headerLine,
                    file = FiltWigOut,
                    sep = "\n",
                    append = T
                    )
                
                write_delim(
                        x = FiltAll,
                        path = FiltWigOut,
                        delim = " ",
                        col_names = F,
                        append = T
                )
                cat("\n\n")
                
                #Clean large variable from memory
                # rm(BedFile)
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
