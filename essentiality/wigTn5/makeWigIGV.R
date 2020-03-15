library(tidyverse)

for (seqLib in list.dirs(recursive = F)) {
        seqRoot <- str_split(seqLib, "/", simplify = T)[2]
        print(paste("Processing Sample: ", seqRoot))
        bedIN <- paste(seqLib, "/alignment/", seqRoot, ".bed", sep = '')
        print(paste("Reading in:", bedIN), quote = F)
        BedFile <- read_tsv(bedIN, col_names = F, trim_ws = T)
        names(BedFile) <- c("Chrom",
                            "ChromStart",
                            "ChromEnd",
                            "ReadID",
                            "MPQ",
                            "Strand")
        print(paste("Total Reads: ", length(BedFile$Chrom)))

        PlusStrand <- BedFile %>% filter(Strand == "+")

        MinusStrand <- BedFile %>% filter(Strand == "-")

        print(paste("Unique Insertions Plus Strand: ", length(unique(
                PlusStrand$ChromStart
        ))))
        print(paste("Unique Insertions Minus Strand: ", length(unique(
                MinusStrand$ChromStart
        ))))
        wigOutPath <- paste(seqLib, "/wig", sep = '')
        dir.create(path = wigOutPath, showWarnings = T)

        PlusWig <- as.data.frame(table(PlusStrand$ChromStart+10))
        PWigOut <- paste(wigOutPath, "/", seqRoot, ".plus.wig", sep = "")
        file.create(PWigOut)
        cat("variableStep chrom=NZ_AKIA01000001.1",file=PWigOut,sep="\n")
        write_delim(x = PlusWig,
                    path = PWigOut,
                    delim = " ",
                    col_names = F,
                    append = T)


        MinusWig <-as.data.frame(table(MinusStrand$ChromEnd-10))
        MinusWig$Freq <- MinusWig$Freq*(-1)

        MWigOut <- paste(wigOutPath, "/", seqRoot, ".minus.wig", sep = "")
        file.create(MWigOut)
        cat("variableStep chrom=NZ_AKIA01000001.1",file=MWigOut,sep="\n")
        write_delim(x = MinusWig,
                    path = MWigOut,
                    delim = " ",
                    col_names = F,
                    append = T)

        ObsAll <- ifelse(BedFile$Strand == "+", BedFile$ChromStart+10, BedFile$ChromEnd-10) # Combined offsets
        AllWig <- as.data.frame(table(ObsAll))
        AllWigOut <- paste(wigOutPath, "/", seqRoot, ".all.wig", sep = "")
        file.create(AllWigOut)
        cat("variableStep chrom=NZ_AKIA01000001.1",file=AllWigOut,sep="\n")
        write_delim(x = AllWig,
                    path = AllWigOut,
                    delim = " ",
                    col_names = F,
                    append = T)
}
