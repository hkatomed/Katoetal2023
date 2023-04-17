#########################################################
## 34_calcHBAexonsMouse2.R
## Calculation of HBA along protein-coding genes (mouse)
## Liftover of nucleosomes: 20170912_mouse_liftOver_1.txt (2)
## R:
setwd(HOMEDIR)
setwd("GSE82127_RAW")
beforelift <- read.table(file = "beforelift_mm9_to_mm10.txt")
lifted <- read.table(file = "lifted_mm9_to_mm10.txt")
unlifted <- read.table(file = "unlifted_mm9_to_mm10.txt")
without_otherChr <- read.table(file = "lifted_mm9_to_mm10_without_otherChr.txt")
otherChr_added <- read.table(file = "lifted_mm9_to_mm10_otherChr_added.txt")
otherChr_ls <- read.table(file = "otherChr_ls.txt")
random_ls <- read.table(file = "random_ls.txt")
sum(beforelift$V1)
sum(lifted$V1)
sum(unlifted$V1)
sum(otherChr_ls$V1)
sum(without_otherChr$V1)
sum(otherChr_added$V1)
sum(random_ls$V1)


