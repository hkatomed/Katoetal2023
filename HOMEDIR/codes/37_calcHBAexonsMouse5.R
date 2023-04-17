#########################################################
## 37_calcHBAexonsMouse5.R
## Calculation of HBA along protein-coding genes (mouse)
## Preparation of unique nucleosomes
## generation of mm9 bed: 20211214_nuCpos2_1.txt (1)	
library(Biostrings)
library(parallel)
setwd(HOMEDIR)
setwd("GSE82127_RAW")
uniqueNuc <- read.table(file = "GSM2183909_unique.map_95pc.txt.gz", header = FALSE)
for(i in 1:22){
	cat(i, ",", sep = "")
	if(i < 20)	uniqueNuc$V1[uniqueNuc$V1 == i] <- paste("chr", i, sep = "")
	if(i == 20)	uniqueNuc$V1[uniqueNuc$V1 == i] <- "chrX"
	if(i == 21)	uniqueNuc$V1[uniqueNuc$V1 == i] <- "chrY"
	if(i == 22)	uniqueNuc$V1[uniqueNuc$V1 == i] <- "chrM"
}
uniqueNuc$V3 <- uniqueNuc$V2
colnames(uniqueNuc) <- c("chr", "start", "end")
uniqueNuc$start <- uniqueNuc$start - 1
class(uniqueNuc$start) <- "integer"
setwd(HOMEDIR)
setwd("20211210_mouse_GRCm38_p6")
dir.create("liftOver_unique")
setwd("liftOver_unique")
write.table(uniqueNuc, file = "uniqueNuc_mm9.bed", 
	sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
rm(uniqueNuc)


