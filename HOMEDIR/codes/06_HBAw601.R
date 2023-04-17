#########################################################
## 06_HBAw601.R
## Calculate HBA along Widom 601 (Figure 3C): 20210716_nuCpos2_1.txt
setwd(HOMEDIR)
setwd("20210706_parameters")
setwd("RData")
load(file = "nature11142_s2_67352_147_freqN4.RData")
load(file = "nature11142_s2_67352_147_tranN4_SMA.RData")
load(file = "merged33_67352_147_freqN4.RData")
load(file = "merged33_67352_147_tranN4_SMA.RData")
library(nuCpos)
setwd("/scripts/nuCpos_R")
source("HBA6.R")
load("sysdata.rda")
SPECIES <- c("scH4S47CL", "scH3H4_33L")
ori601 <- "CGGGATCCTAATGACCAAGGAAAGCATGATTCTTCACACCGAGTTCATCCCTTATGTGATGGACCCTATACGCGGCCGCCCTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGTGCATGTATTGAACAGCGACCTTGCCGGTGCCAGTCGGATAGTGTTCCGAGCTCCC"
HBAdata <- data.frame(pos = 74:209, stringsAsFactors = FALSE)
for(s in 1:length(SPECIES)){
	species <- SPECIES[s]
	HBAdata[[species]] <- numeric(length = 136)
	for(i in 1:136){
		seq <- substr(x = ori601, start = i, stop = i+146)
		HBAdata[[species]][i] <- HBA6(seq, species = species, silent = TRUE)
	}
}
HBAdata2 <- HBAdata
for(s in 1:length(SPECIES)){
	species <- SPECIES[s]
	HBAdata2[[species]] <- HBAdata2[[species]] - mean(HBAdata2[[species]])
}
HBAdata3 <- HBAdata
for(s in 1:length(SPECIES)){
	species <- SPECIES[s]
	HBAdata3[[species]] <- scale(HBAdata3[[species]])
}
ymax <- max(HBAdata3[, c(2:8)])
ymin <- min(HBAdata3[, c(2:8)])
ylim <- c(ymin, ymax)
drawHBA2 <- function(species){
	main <- paste("scH3H4_33L vs ", species, sep = "")
	plot(x = 74:209, y = HBAdata3$scH3H4_33L, type = "n", lwd = 2, 
		xlab = "Dyad position (bp)", ylab = "Standardized HBA", ylim = ylim,
		main = main)
	abline(v = 154, lwd = 2, col = "orange")
	lines(x = 74:209, y = HBAdata3[[species]], type = "l", lwd = 2, col = rgb(0,0,1, alpha = 0.5), lty = 1)
	lines(x = 74:209, y = HBAdata3$scH3H4_33L, type = "l", lwd = 2, col = rgb(1,0,0, alpha = 0.5), lty = 1)
	legend(x = 160, y = -1.3, legend = "scH3H4_33L", bty = "n", lwd = 2, col = rgb(1,0,0, alpha = 0.5), lty = 1)
	legend(x = 160, y = -1.8, legend = species, bty = "n", lwd = 2, col = rgb(0,0,1, alpha = 0.5), lty = 1)
	box(lwd = 2)
	axis(side = 1, lwd.ticks = 2)
	axis(side = 2, lwd.ticks = 2)
}
drawHBA2("scH4S47CL")


