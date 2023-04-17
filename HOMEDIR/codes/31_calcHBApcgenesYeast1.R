#########################################################
## 31_calcHBApcgenesYeast1.R
## Calculation of HBA along protein-coding genes (yeast)
## Calculation: 20210906_nuCpos2_1.txt
setwd(HOMEDIR)
setwd("20210706_parameters/scripts/nuCpos_src")
dyn.load("LRA_1.so")
setwd(HOMEDIR)
setwd("20210706_parameters/scripts/nuCpos_R")
source("LRA1.R")
load("sysdata.rda")
setwd(HOMEDIR)
setwd("20210706_parameters")
setwd("RData")
load(file = "nature11142_s2_147_freqN4R.RData")
load(file = "merged33_67352_147_freqN4R.RData")
load(file = "nature11142_s2_147_freqN4R2.RData")
load(file = "merged33_67352_147_freqN4R2.RData")	
load(file = "nature11142_s2_147_freqN4R3.RData")
load(file = "merged33_67352_147_freqN4R3.RData")	
load(file = "nature11142_s2_147_freqN4R4.RData")
load(file = "merged33_67352_147_freqN4R4.RData")	
load(file = "nature11142_s2_147_freqN4R5.RData")
load(file = "merged33_67352_147_freqN4R5.RData")	
load(file = "merged33_67352_147_freqN4SA.RData")	
load(file = "merged33_67352_147_freqN4SC.RData")
load(file = "merged33_67352_147_tranN4_SMA.RData")
setwd(HOMEDIR)
setwd("20210706_parameters")
setwd("RData")
load(file = "sc_genome_2011.RData")
setwd(HOMEDIR)
setwd("20210708_Chereji_HBA")
setwd("RData")
load(file = "Chereji_over1500.RData")
get.LRA1 <- function(i, half = "L"){
	strand <- Chereji.over1500$Strand[i]
	chr <- Chereji.over1500$Chr[i]
	if(strand == 1){
		ORF.left <- Chereji.over1500$plus1[i] -500-73
		ORF.right <- Chereji.over1500$plus1[i] +1000+73		
		seq1646 <- sc.genome.2011[[chr]][ORF.left:ORF.right]
	}
	if(strand == -1){
		ORF.left <- Chereji.over1500$plus1[i] -1000-73		
		ORF.right <- Chereji.over1500$plus1[i] +500+73
		seq1646 <- reverseComplement(sc.genome.2011[[chr]][ORF.left:ORF.right])
	}
	rowdata <- numeric(length = 1501)
	for(s in 1:1501){
		temp <- LRA1(seq1646[s:(s+146)], species = species, silent = TRUE)
		if(half == "L")	rowdata[s] <- temp[1]
		if(half == "R")	rowdata[s] <- temp[2]
	}
	return(rowdata)
}
setwd(HOMEDIR)
dir.create("20210906_Chereji_HBA")
setwd("20210906_Chereji_HBA")
dir.create("RData")
setwd("RData")
species <- "scH3H4_33L"
matrix1646.LRA1.L <- matrix(nrow = 1501, ncol = nrow(Chereji.over1500))	
colnames(matrix1646.LRA1.L) <- Chereji.over1500$ORF
matrix1646.LRA1.L[,1:nrow(Chereji.over1500)] <- 
	unlist(mcMap(f = get.LRA1, i = 1:nrow(Chereji.over1500), half = "L", mc.cores = 36), use.names = FALSE)
matrix1646.LRA1.R <- matrix(nrow = 1501, ncol = nrow(Chereji.over1500))	
colnames(matrix1646.LRA1.R) <- Chereji.over1500$ORF
matrix1646.LRA1.R[,1:nrow(Chereji.over1500)] <- 
	unlist(mcMap(f = get.LRA1, i = 1:nrow(Chereji.over1500), half = "R", mc.cores = 36), use.names = FALSE)
setwd(HOMEDIR)
setwd("20210906_Chereji_HBA")
setwd("RData")
save(matrix1646.LRA1.L, file = "matrix1646.LRA1.L_Chereji.RData")
save(matrix1646.LRA1.R, file = "matrix1646.LRA1.R_Chereji.RData")
library(nuCpos)
setwd(HOMEDIR)
setwd("20210706_parameters")
setwd("RData")
load(file = "merged33_67352_147_tranN4_SMA.RData")
load(file = "merged33_67352_147_freqN4.RData")
setwd(HOMEDIR)
setwd("20210706_parameters/scripts/nuCpos_R")
source("HBA5.R")
load("sysdata.rda")
get.HBA <- function(i){
	strand <- Chereji.over1500$Strand[i]
	chr <- Chereji.over1500$Chr[i]
	if(strand == 1){
		ORF.left <- Chereji.over1500$plus1[i] -500-73
		ORF.right <- Chereji.over1500$plus1[i] +1000+73
		seq1646 <- sc.genome.2011[[chr]][ORF.left:ORF.right]
	}
	if(strand == -1){
		ORF.left <- Chereji.over1500$plus1[i] -1000-73
		ORF.right <- Chereji.over1500$plus1[i] +500+73
		seq1646 <- reverseComplement(sc.genome.2011[[chr]][ORF.left:ORF.right])
	}
	rowdata <- numeric(length = 1501)
	for(s in 1:1501){
		rowdata[s] <- HBA5(seq1646[s:(s+146)], species = species, silent = TRUE)
	}
	return(rowdata)
}
species <- "scH3H4_33L"
matrix1646 <- matrix(nrow = 1501, ncol = nrow(Chereji.over1500))
colnames(matrix1646) <- Chereji.over1500$ORF
matrix1646[,1:nrow(Chereji.over1500)] <- 
	unlist(mcMap(f = get.HBA, i = 1:nrow(Chereji.over1500), mc.cores = 36), use.names = FALSE)
setwd(HOMEDIR)
setwd("20210906_Chereji_HBA")
setwd("RData")
save(matrix1646, file = "matrix1646_Chereji.RData")


