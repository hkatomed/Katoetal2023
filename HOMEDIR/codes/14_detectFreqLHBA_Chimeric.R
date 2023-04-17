#########################################################
## 14_detectFreqLHBA_Chimeric.R
## Analyze the relationship between detection frequency and local HBA
## with chimeric models (Figure 2): 20210721_nuCpos2_7.txt
setwd(HOMEDIR)
setwd("20210706_parameters")
setwd("RData")
load(file = "merged33_67352_147_freqN4.RData")
load(file = "merged33_67352_147_tranN4_SMA.RData")
load(file = "merged33_67352_147_freqN4SA.RData")
load(file = "merged33_67352_147_freqN4SB.RData")
load(file = "merged33_67352_147_freqN4SC.RData")
load(file = "merged33_67352_147_freqN4SD.RData")
load(file = "merged33_67352_147_freqN4SE.RData")
load(file = "merged33_67352_147_freqN4SF.RData")
load(file = "merged33_67352_147_freqN4SG.RData")
load(file = "merged33_67352_147_freqN4SH.RData")
load(file = "merged33_67352_147_freqN4SI.RData")
load(file = "merged33_67352_147_freqN4SJ.RData")
load(file = "merged33_67352_147_freqN4SK.RData")
load(file = "merged33_67352_147_freqN4SL.RData")
load(file = "merged33_67352_147_freqN4SM.RData")
setwd(HOMEDIR)
setwd("20210706_nuCpos2")
setwd("scripts/nuCpos_R")
source("HBA6.R")
source("localHBA2.R")
load("sysdata.rda")
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
load(file = "Chereji_w107_add3.RData")
load(file = "Brogaard_Uniq_add3.RData")
library(Biostrings)
library(nuCpos)
library(parallel)
Chereji_w107$HBA_sc33L <- unlist(mcMap(f = HBA6, inseq = Chereji_w107$seq,
					species = "scH3H4_33L", silent = TRUE, mc.cores = 24), 
					use.names = FALSE)
Brogaard_Uniq$HBA_sc33L <- unlist(mcMap(f = HBA6, inseq = Brogaard_Uniq$seq,
					species = "scH3H4_33L", silent = TRUE, mc.cores = 24), 
					use.names = FALSE)
Chereji_w107$lA_sc33L <- as.numeric(NA)
Chereji_w107$lB_sc33L <- as.numeric(NA)
Chereji_w107$lC_sc33L <- as.numeric(NA)
Chereji_w107$lD_sc33L <- as.numeric(NA)
Chereji_w107$lE_sc33L <- as.numeric(NA)
Chereji_w107$lF_sc33L <- as.numeric(NA)
Chereji_w107$lG_sc33L <- as.numeric(NA)
Chereji_w107$lH_sc33L <- as.numeric(NA)
Chereji_w107$lI_sc33L <- as.numeric(NA)
Chereji_w107$lJ_sc33L <- as.numeric(NA)
Chereji_w107$lK_sc33L <- as.numeric(NA)
Chereji_w107$lL_sc33L <- as.numeric(NA)
Chereji_w107$lM_sc33L <- as.numeric(NA)
Chereji_w107[, 36:48] <- matrix(unlist(
			mcMap(f = localHBA2, inseq = Chereji_w107$seq,
				species = "scH3H4_33L", silent = TRUE, mc.cores = 24), 
				use.names = FALSE), byrow = TRUE, ncol = 13)
Brogaard_Uniq$lA_sc33L <- as.numeric(NA)
Brogaard_Uniq$lB_sc33L <- as.numeric(NA)
Brogaard_Uniq$lC_sc33L <- as.numeric(NA)
Brogaard_Uniq$lD_sc33L <- as.numeric(NA)
Brogaard_Uniq$lE_sc33L <- as.numeric(NA)
Brogaard_Uniq$lF_sc33L <- as.numeric(NA)
Brogaard_Uniq$lG_sc33L <- as.numeric(NA)
Brogaard_Uniq$lH_sc33L <- as.numeric(NA)
Brogaard_Uniq$lI_sc33L <- as.numeric(NA)
Brogaard_Uniq$lJ_sc33L <- as.numeric(NA)
Brogaard_Uniq$lK_sc33L <- as.numeric(NA)
Brogaard_Uniq$lL_sc33L <- as.numeric(NA)
Brogaard_Uniq$lM_sc33L <- as.numeric(NA)
Brogaard_Uniq[, 36:48] <- matrix(unlist(
			mcMap(f = localHBA2, inseq = Brogaard_Uniq$seq,
				species = "scH3H4_33L", silent = TRUE, mc.cores = 24), 
				use.names = FALSE), byrow = TRUE, ncol = 13)
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
save(Chereji_w107, file = "Chereji_w107_add4.RData")
save(Brogaard_Uniq, file = "Brogaard_Uniq_add4.RData")
temp <- Chereji_w107
temp <- subset(temp, dyadScoreGroup != "rmHigh" & dyadScoreGroup != "rmLow")
meanHBA <- numeric(length = 9)
meanLocalHBA <- data.frame(matrix(data = 0, nrow = 13, ncol = 9), stringsAsFactors = FALSE)
for(g in 1:9){
	sg <- paste("sg", g, sep = "")
	temp2 <- subset(temp, dyadScoreGroup == sg)
	names(meanHBA)[g] <- sg
	colnames(meanLocalHBA)[g] <- sg
	meanHBA[g] <- mean(temp2$HBA_sc33L, na.rm = TRUE)
	meanLocalHBA[,g] <- colMeans(temp2[,36:48], na.rm = TRUE)
	if(g == 1) rownames(meanLocalHBA) <- names(colMeans(temp2[,36:48], na.rm = TRUE))
}
Chereji.meanHBA <- meanHBA
Chereji.meanLocalHBA <- meanLocalHBA
temp <- Brogaard_Uniq
temp <- subset(temp, dyadScoreGroup != "rmHigh" & dyadScoreGroup != "rmLow")
meanHBA <- numeric(length = 9)
meanLocalHBA <- data.frame(matrix(data = 0, nrow = 13, ncol = 9), stringsAsFactors = FALSE)
for(g in 1:9){
	sg <- paste("sg", g, sep = "")
	temp2 <- subset(temp, dyadScoreGroup == sg)
	names(meanHBA)[g] <- sg
	colnames(meanLocalHBA)[g] <- sg
	meanHBA[g] <- mean(temp2$HBA_sc33L, na.rm = TRUE)
	meanLocalHBA[,g] <- colMeans(temp2[,36:48], na.rm = TRUE)
	if(g == 1) rownames(meanLocalHBA) <- names(colMeans(temp2[,36:48], na.rm = TRUE))
}
Brogaard.meanHBA <- meanHBA
Brogaard.meanLocalHBA <- meanLocalHBA
ylim <- c(min(Chereji.meanHBA, Brogaard.meanHBA), max(Chereji.meanHBA, Brogaard.meanHBA))
plot(x = 1:9, y = numeric(length = 9), ylim = ylim, ylab = "Mean HBA (scH3H4_33L)", 
		xaxt = "n", yaxt = "n", type = "n")
lines(x = 1:9, y = Chereji.meanHBA, col = "red", lwd = 3)
lines(x = 1:9, y = Brogaard.meanHBA, col = "blue", lwd = 3)
axis(1, at = 1:9, lwd = 3)
axis(2, lwd = 3)
box(lwd = 3)
legend(x = 6, y = 2.9, legend = "Brogaard", lwd = 3, col = "blue", bty = "n")
legend(x = 6, y = 2.75, legend = "Chereji", lwd = 3, col = "red", bty = "n")
library(RColorBrewer)
COL <- brewer.pal(9, "Spectral")
ylim <- c(min(Chereji.meanLocalHBA, Brogaard.meanLocalHBA), 
	max(Chereji.meanLocalHBA, Brogaard.meanLocalHBA))
main <- "Chereji local HBA (scH3H4_33L)"
plot(x = 1:13, y = numeric(length = 13), ylim = ylim, ylab = "Mean local HBA (scH3H4_33L)", 
		xaxt = "n", yaxt = "n", type = "n", main = main)
axis(1, at = 1:13, labels = LETTERS[1:13], lwd = 3)
axis(2, lwd = 3)
for(g in 1:9){
	lines(x = 1:13, y = Chereji.meanLocalHBA[,g], lwd = 3, col = COL[g])
	legend(x = 11.2, y = 0.41 - 0.02*g, legend = paste("sg", g, sep = ""), 
		lwd = 3, col = COL[g], bty = "n")
}
box(lwd = 3)
main <- "Brogaard local HBA (scH3H4_33L)"
plot(x = 1:13, y = numeric(length = 13), ylim = ylim, ylab = "Mean local HBA (scH3H4_33L)", 
		xaxt = "n", yaxt = "n", type = "n", main = main)
axis(1, at = 1:13, labels = LETTERS[1:13], lwd = 3)
axis(2, lwd = 3)
for(g in 1:9){
	lines(x = 1:13, y = Brogaard.meanLocalHBA[,g], lwd = 3, col = COL[g])
	legend(x = 11.2, y = 0.205 - 0.02*g, legend = paste("sg", g, sep = ""), 
		lwd = 3, col = COL[g], bty = "n")
}
box(lwd = 3)
ylim <- c(-0.5801816, 1.6358595)
main <- "Chereji local HBA (scH3H4_33L)"
plot(x = 1:13, y = numeric(length = 13), ylim = ylim, ylab = "Mean local HBA (scH3H4_33L)", 
		xaxt = "n", yaxt = "n", type = "n", main = main)
axis(1, at = 1:13, labels = LETTERS[1:13], lwd = 3)
axis(2, lwd = 3)
for(g in 1:9){
	lines(x = 1:13, y = Chereji.meanLocalHBA[,g], lwd = 3, col = COL[g])
	legend(x = 11.2, y = 1.8 - 0.12*g, legend = paste("sg", g, sep = ""), 
		lwd = 3, col = COL[g], bty = "n")
}
box(lwd = 3)
main <- "Brogaard local HBA (scH3H4_33L)"
plot(x = 1:13, y = numeric(length = 13), ylim = ylim, ylab = "Mean local HBA (scH3H4_33L)", 
		xaxt = "n", yaxt = "n", type = "n", main = main)
axis(1, at = 1:13, labels = LETTERS[1:13], lwd = 3)
axis(2, lwd = 3)
for(g in 1:9){
	lines(x = 1:13, y = Brogaard.meanLocalHBA[,g], lwd = 3, col = COL[g])
	legend(x = 11.2, y = 1.8 - 0.12*g, legend = paste("sg", g, sep = ""), 
		lwd = 3, col = COL[g], bty = "n")
}
box(lwd = 3)


