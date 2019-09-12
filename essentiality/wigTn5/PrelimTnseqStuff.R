library(tidyverse)
library(vegan)

BedFile <- read_tsv(file = "MESO.bed", col_names = F, trim_ws = T)
names(BedFile) <- c("Chrom", "ChromStart", "ChromEnd", "ReadID", "MPQ", "Strand")

#Split by Strand
PlusStrand <- BedFile %>% filter(Strand == "+")
MinusStrand <- BedFile %>% filter(Strand == "-")

#Sanity Checks
length(unique(PlusStrand$ChromStart))
length(PlusStrand$ChromStart)
PlusWig <- as.data.frame(table(PlusStrand$ChromStart))
write_delim(x = PlusWig, path = "W1.Plus.wig", delim = " ")
length(unique(MinusStrand$ChromEnd))
length(MinusStrand$ChromEnd)
MinusWig <-as.data.frame(table(MinusStrand$ChromEnd))
MinusWig$Freq <- MinusWig$Freq*(-1)
write_delim(x = MinusWig, path = "W1.Minus.wig", delim = " ")

#Total Unique Insertions: 43942 + 44834 = 88,776 

length(BedFile$Chrom)
#Total Reads: 340,492
length(ifelse(BedFile$Strand == "+", BedFile$ChromStart, BedFile$ChromEnd))
#Length: 340492
length(unique(ifelse(BedFile$Strand == "+", BedFile$ChromStart, BedFile$ChromEnd)))
#Length: 88558
#Meaning that 88776 uniques when split by strand, 88558 when not split by strand,
#therefore 218 insertions in Fwd and Rev share the same insertion site.
ObsAll <- ifelse(BedFile$Strand == "+", BedFile$ChromStart, BedFile$ChromEnd)
CountTable <- as.data.frame(table(ifelse(BedFile$Strand == "+", BedFile$ChromStart, BedFile$ChromEnd)))
sum(CountTable$Freq)
#Total Counts: 340492
names(CountTable) <- c("Pos", "Count")
CountTable$Pos <- paste("Pos", CountTable$Pos, sep = '')
write_tsv(x = CountTable, path = "4815.1.counts.tsv")

hist(x = CountTable$Count, breaks = 1000, main = "Mutant Abundance Distrbution", xlab = "Mutant Abundance", ylab = "Unique Mutants", xlim = c(1,50)
length(table(CountTable$Count))


