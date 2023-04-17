#########################################################
## 07_localHBAw601.R
## Calculate local HBA along Widom 601 (Figure 3D): 20210716_nuCpos2_2.txt
setwd(HOMEDIR)
setwd("20210706_parameters")
setwd("scripts/nuCpos_R")
source("localHBA2.R")
load("sysdata.rda")
setwd("../RData")
load(file = "nature11142_s2_67352_147_freqN4SA.RData")
load(file = "nature11142_s2_67352_147_freqN4SB.RData")
load(file = "nature11142_s2_67352_147_freqN4SC.RData")
load(file = "nature11142_s2_67352_147_freqN4SD.RData")
load(file = "nature11142_s2_67352_147_freqN4SE.RData")
load(file = "nature11142_s2_67352_147_freqN4SF.RData")
load(file = "nature11142_s2_67352_147_freqN4SG.RData")
load(file = "merged33_67352_147_freqN4SA.RData")
load(file = "merged33_67352_147_freqN4SB.RData")
load(file = "merged33_67352_147_freqN4SC.RData")
load(file = "merged33_67352_147_freqN4SD.RData")
load(file = "merged33_67352_147_freqN4SE.RData")
load(file = "merged33_67352_147_freqN4SF.RData")
load(file = "merged33_67352_147_freqN4SG.RData")
load(file = "nature11142_s2_67352_147_freqN4SH.RData")
load(file = "nature11142_s2_67352_147_freqN4SI.RData")
load(file = "nature11142_s2_67352_147_freqN4SJ.RData")
load(file = "nature11142_s2_67352_147_freqN4SK.RData")
load(file = "nature11142_s2_67352_147_freqN4SL.RData")
load(file = "nature11142_s2_67352_147_freqN4SM.RData")
load(file = "merged33_67352_147_freqN4SH.RData")
load(file = "merged33_67352_147_freqN4SI.RData")
load(file = "merged33_67352_147_freqN4SJ.RData")
load(file = "merged33_67352_147_freqN4SK.RData")
load(file = "merged33_67352_147_freqN4SL.RData")
load(file = "merged33_67352_147_freqN4SM.RData")
load(file = "nature11142_s2_67352_147_freqN4.RData")
load(file = "nature11142_s2_67352_147_tranN4_SMA.RData")
load(file = "merged33_67352_147_freqN4.RData")
load(file = "merged33_67352_147_tranN4_SMA.RData")
SPECIES <- c("scH4S47CL", "scH3H4_33L")
ori601 <- "CGGGATCCTAATGACCAAGGAAAGCATGATTCTTCACACCGAGTTCATCCCTTATGTGATGGACCCTATACGCGGCCGCCCTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGTGCATGTATTGAACAGCGACCTTGCCGGTGCCAGTCGGATAGTGTTCCGAGCTCCC"
lHBAdata <- list()
for(s in 1:length(SPECIES)){
	species <- SPECIES[s]
	rawdata <- data.frame(pos = 74:209, lHBA_A = as.numeric(NA), lHBA_B = as.numeric(NA), 
				lHBA_C = as.numeric(NA), lHBA_D = as.numeric(NA), 
				lHBA_E = as.numeric(NA), lHBA_F = as.numeric(NA), 
				lHBA_G = as.numeric(NA), lHBA_H = as.numeric(NA), 
				lHBA_I = as.numeric(NA), lHBA_J = as.numeric(NA), 
				lHBA_K = as.numeric(NA), lHBA_L = as.numeric(NA), 
				lHBA_M = as.numeric(NA), stringsAsFactors = FALSE)
	for(i in 1:136){
		seq <- substr(x = ori601, start = i, stop = i+146)
		rawdata[i,2:14] <- localHBA2(seq, species = species, silent = TRUE)
	}
	lHBAdata[[species]] <- rawdata
}
setwd(HOMEDIR)
dir.create("20210716_Widom_localHBA")
setwd("20210716_Widom_localHBA")
library(rasterVis)
COL = BuRdTheme()$regions$col
for(s in 1:length(SPECIES)){
	species <- SPECIES[s]
	grid <- as.matrix(lHBAdata[[species]][, 2:14])
	grid <- grid[, ncol(grid):1]
	cat("# ", species, ": ", max(grid), ", ", min(grid), "\n")
	gridNum <- 6
	gridNew <- matrix(nrow = nrow(grid), ncol = ncol(grid)*gridNum)
	for(i in 1:ncol(grid)) {
		row.range <- ((i-1)*gridNum + 1):((i-1)*gridNum + gridNum)
		gridNew[, row.range] <- grid[,i]
	}	
	out <- levelplot(x = gridNew, useRaster = TRUE, at = seq(-9, 5, 1),
		colorkey = list(at = seq(-9, 5, 0.5), col = COL), col.regions = COL, 
		scales = list(x = list(at = seq(7, 127, 20), labels = seq(80, 200, 20)), 
				y = list(at = seq(gridNum/2, nrow(gridNew), gridNum), 
				labels = LETTERS[13:1])), 
		xlab = "Dyad position (bp)", 
		ylab = "Segments")
	filename <- paste("Widom601_282_localHBA_", species, ".jpeg", sep = "")
	jpeg(file = filename, width = 900, height = 600, quality = 200, res = 200)
	print(out)
	dev.off()
	filename <- paste("Widom601_282_localHBA_", species, ".tiff", sep = "")
	tiff(file = filename, width = 900, height = 600, res = 200)
	print(out)
	dev.off()
}


