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

        print(paste("Total Reads: ", length(BedFile$Chrom)))
        cat("\n")
        PlusStrand <- BedFile %>% filter(Strand == "+")
        MinusStrand <- BedFile %>% filter(Strand == "-")

        wigOutPath <- paste(seqLib, "/wig", sep = '')
        dir.create(path = wigOutPath, showWarnings = T)


        #Plus Strand Wig Stats and Wig File
        PlusWig <- as.data.frame(table(PlusStrand$ChromStart+10))
        PWigOut <- paste(wigOutPath, "/", seqRoot, ".plus.wig", sep = "")
        print(paste("Printing .wig file:", PWigOut))
        file.create(PWigOut)
        print(paste("Unique Insertions Plus Strand: ", length(
                PlusWig$Freq
        )))
        cat("Summary Statistics:\n")
        print(describe(PlusWig$Freq, trim = 0.05, IQR = T))
        PlusWig$Var1 <- as.integer(levels(PlusWig$Var1))
        cat("variableStep chrom=CP051772",file=PWigOut,sep="\n")
        write_delim(x = PlusWig,
                    path = PWigOut,
                    delim = " ",
                    col_names = F,
                    append = T)
        cat("\n\n")


        #Minus Strand Wig Stats and Wig File
        MinusWig <-as.data.frame(table(MinusStrand$ChromEnd-10))
        MinusWigforDescribe <- MinusWig
        MinusWig$Freq <- MinusWig$Freq*(-1)
        MWigOut <- paste(wigOutPath, "/", seqRoot, ".minus.wig", sep = "")
        print(paste("Printing .wig file:", MWigOut))
        file.create(MWigOut)
        print(paste("Unique Insertions Minus Strand: ", length(
                MinusWigforDescribe$Freq
        )))
        cat("Summary Statistics:\n")
        print(describe(MinusWigforDescribe$Freq, trim = 0.05, IQR = T))
        MinusWig$Var1 <- as.integer(levels(MinusWig$Var1))
        cat("variableStep chrom=CP051772",file=MWigOut,sep="\n")
        write_delim(x = MinusWig,
                    path = MWigOut,
                    delim = " ",
                    col_names = F,
                    append = T)
        cat("\n\n")


        ObsAll <- ifelse(BedFile$Strand == "+", BedFile$ChromStart+10, BedFile$ChromEnd-10) # Combined offsets
        AllWig <- as.data.frame(table(ObsAll))
        AllWigOut <- paste(wigOutPath, "/", seqRoot, ".all.wig", sep = "")
        print(paste("Printing .wig file:", AllWigOut))
        file.create(AllWigOut)
        print(paste("Unique Insertions In Total: ", length(
                AllWig$Freq
        )))
        cat("Summary Statistics:\n")
        print(describe(AllWig$Freq, trim = 0.05, IQR = T))
        AllWig$ObsAll <- as.integer(levels(AllWig$ObsAll))
        cat("variableStep chrom=CP051772",file=AllWigOut,sep="\n")
        write_delim(x = AllWig,
                    path = AllWigOut,
                    delim = " ",
                    col_names = F,
                    append = T)
        cat("\n\n")


        FiltAll <- AllWig %>% filter(Freq>1) #remove singletons
        FiltWigOut <- paste(wigOutPath, "/", seqRoot, ".filt.1.wig", sep = "")
        print(paste("Printing .wig file:", FiltWigOut))
        file.create(FiltWigOut)
        print(paste("Unique Insertions In Total: ", length(
                FiltAll$Freq
        )))
        cat("Summary Statistics:\n")
        print(describe(FiltAll$Freq, trim = 0.05, IQR = T))
        cat("variableStep chrom=CP051772",file=FiltWigOut,sep="\n")
        write_delim(x = FiltAll,
                    path = FiltWigOut,
                    delim = " ",
                    col_names = F,
                    append = T)
        cat("\n\n")


        FiltAll2 <- AllWig %>% filter(Freq>2) #remove singletons and doubletons
        FiltWigOut2 <- paste(wigOutPath, "/", seqRoot, ".filt.2.wig", sep = "")
        print(paste("Printing .wig file:", FiltWigOut2))
        file.create(FiltWigOut2)
        print(paste("Unique Insertions In Total: ", length(
                FiltAll2$Freq
        )))
        cat("Summary Statistics:\n")
        print(describe(FiltAll2$Freq, trim = 0.05, IQR = T))
        cat("variableStep chrom=CP051772",file=FiltWigOut2,sep="\n")
        write_delim(x = FiltAll2,
                    path = FiltWigOut2,
                    delim = " ",
                    col_names = F,
                    append = T)
        cat("\n\n")

        #Clean large variable from memory
        rm(BedFile)
        rm(PlusStrand)
        rm(MinusStrand)
        rm(PlusWig)
        rm(MinusWig)
        rm(ObsAll)
        rm(AllWig)
        rm(FiltAll)
        rm(FiltAll2)

}
