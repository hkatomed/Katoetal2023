#########################################################
## 43_calcHBAexonsMouse11.R
## Calculation of HBA along protein-coding genes (mouse)
## Calculation: 20211217_nuCpos2_4.txt
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
library(nuCpos)
setwd(HOMEDIR)
setwd("20210706_parameters")
setwd("RData")
load(file = "merged33_67352_147_tranN4_SMA.RData")
load(file = "merged33_67352_147_freqN4.RData")
setwd(HOMEDIR)
setwd("20210706_parameters/scripts/nuCpos_R")
source("HBA5.R")
setwd(HOMEDIR)
setwd("20211210_mouse_GRCm38_p6")
setwd("RData")
load(file = "nucSeq_over500_CentralNuc.RData")
setwd(HOMEDIR)
dir.create("20211215_mm10_HBA")
setwd("20211215_mm10_HBA")
dir.create("RData")
setwd("RData")
HBAlist_CentralNuc <- list()
HBAlist_CentralNuc[["chimeric"]] <- list()
HBAlist_CentralNuc[["sc"]] <- list()
HBAmatrix <- matrix(data = as.numeric(NA), nrow = length(nucSeq_over500_CentralNuc), ncol = 1001)
HBAlist_CentralNuc[["chimeric"]][["HBA"]] <- HBAmatrix
HBAlist_CentralNuc[["chimeric"]][["L"]] <- HBAmatrix
HBAlist_CentralNuc[["chimeric"]][["R"]] <- HBAmatrix
HBAlist_CentralNuc[["sc"]][["HBA"]] <- HBAmatrix
HBAlist_CentralNuc[["sc"]][["L"]] <- HBAmatrix
HBAlist_CentralNuc[["sc"]][["R"]] <- HBAmatrix
for(i in 1:length(nucSeq_over500_CentralNuc)){
	if(i%%20 == 0) cat(i, ",", sep = "")
	SEQ <- nucSeq_over500_CentralNuc[[i]]
	testSEQs <- character(length = 1001)
	for(s in 1:1001) testSEQs[s] <- as.character(SEQ[s:(s+146)])
	temp <- unlist(mcMap(f = HBA5, inseq = testSEQs, 
				species = "scH3H4_33L", silent = TRUE, 
				mc.cores = 24), use.names = FALSE)
	HBAlist_CentralNuc[["chimeric"]][["HBA"]][i,] <- temp
	temp <- mcMap(f = LRA1, inseq = testSEQs, 
				species = "scH3H4_33L", silent = TRUE, 
				mc.cores = 24)
	tempL <- numeric(length = 1001)
	tempR <- numeric(length = 1001)
	for(t in 1:1001){
		tempL[t] <- temp[[t]][1]
		tempR[t] <- temp[[t]][2]
	}
	HBAlist_CentralNuc[["chimeric"]][["L"]][i,] <- tempL
	HBAlist_CentralNuc[["chimeric"]][["R"]][i,] <- tempR
	temp <- unlist(mcMap(f = HBA5, inseq = testSEQs, 
				species = "sc", silent = TRUE, 
				mc.cores = 24), use.names = FALSE)
	HBAlist_CentralNuc[["sc"]][["HBA"]][i,] <- temp

	temp <- mcMap(f = LRA1, inseq = testSEQs, 
				species = "sc", silent = TRUE, 
				mc.cores = 24)
	tempL <- numeric(length = 1001)
	tempR <- numeric(length = 1001)
	for(t in 1:1001){
		tempL[t] <- temp[[t]][1]
		tempR[t] <- temp[[t]][2]
	}
	HBAlist_CentralNuc[["sc"]][["L"]][i,] <- tempL
	HBAlist_CentralNuc[["sc"]][["R"]][i,] <- tempR
}
setwd(HOMEDIR)
setwd("20211210_mouse_GRCm38_p6")
setwd("RData")
save(HBAlist_CentralNuc, file = "nucSeq_over500_CentralNuc_HBAlist_CentralNuc.RData")


