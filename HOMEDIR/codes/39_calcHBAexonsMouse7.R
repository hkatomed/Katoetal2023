#########################################################
## 39_calcHBAexonsMouse7.R
## Calculation of HBA along protein-coding genes (mouse)
## Preparation of unique nucleosomes
## Prepare 1-based nucleosome data: 20211214_nuCpos2_1.txt (3)		uniqueNuc_mm10.RData
setwd(HOMEDIR)
setwd("20211210_mouse_GRCm38_p6")
setwd("liftOver_unique")
uniqueNuc_mm10 <- read.table(file = "uniqueNuc_mm10.bed")
colnames(uniqueNuc_mm10)[1] <- "chr"
colnames(uniqueNuc_mm10)[3] <- "pos"
uniqueNuc_mm10$V2 <- NULL
setwd("../RData")
save(uniqueNuc_mm10, file = "uniqueNuc_mm10.RData")


