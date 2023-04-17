#########################################################
## 12_detectFreqLHBA_H4S47C.R
## Analyze the relationship between detection frequency and local HBA
## with H4-S47C–based models (Figure 2): 20210721_nuCpos2_5.txt
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
load(file = "Chereji_w107_add2.RData")
load(file = "Brogaard_Uniq_add2.RData")
temp <- Chereji_w107
temp <- subset(temp, dyadScoreGroup != "rmHigh" & dyadScoreGroup != "rmLow")
meanHBA <- numeric(length = 9)
meanLocalHBA <- data.frame(matrix(data = 0, nrow = 13, ncol = 9), stringsAsFactors = FALSE)
for(g in 1:9){
	sg <- paste("sg", g, sep = "")
	temp2 <- subset(temp, dyadScoreGroup == sg)
	names(meanHBA)[g] <- sg
	colnames(meanLocalHBA)[g] <- sg
	meanHBA[g] <- mean(temp2$HBA_sc, na.rm = TRUE)
	meanLocalHBA[,g] <- colMeans(temp2[,6:18], na.rm = TRUE)
	if(g == 1) rownames(meanLocalHBA) <- names(colMeans(temp2[,6:18], na.rm = TRUE))
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
	meanHBA[g] <- mean(temp2$HBA_sc, na.rm = TRUE)
	meanLocalHBA[,g] <- colMeans(temp2[,6:18], na.rm = TRUE)
	if(g == 1) rownames(meanLocalHBA) <- names(colMeans(temp2[,6:18], na.rm = TRUE))
}
Brogaard.meanHBA <- meanHBA
Brogaard.meanLocalHBA <- meanLocalHBA
ylim <- c(min(Chereji.meanHBA, Brogaard.meanHBA), max(Chereji.meanHBA, Brogaard.meanHBA))
plot(x = 1:9, y = numeric(length = 9), ylim = ylim, ylab = "Mean HBA (sc)", 
		xaxt = "n", yaxt = "n", type = "n")
lines(x = 1:9, y = Chereji.meanHBA, col = "red", lwd = 3)
lines(x = 1:9, y = Brogaard.meanHBA, col = "blue", lwd = 3)
axis(1, at = 1:9, lwd = 3)
axis(2, lwd = 3)
box(lwd = 3)
legend(x = 6, y = 4.3, legend = "Brogaard", lwd = 3, col = "blue", bty = "n")
legend(x = 6, y = 3.9, legend = "Chereji", lwd = 3, col = "red", bty = "n")
library(RColorBrewer)
COL <- brewer.pal(9, "Spectral")
ylim <- c(min(Chereji.meanLocalHBA, Brogaard.meanLocalHBA), 
	max(Chereji.meanLocalHBA, Brogaard.meanLocalHBA))
main <- "Chereji local HBA (sc)"
plot(x = 1:13, y = numeric(length = 13), ylim = ylim, ylab = "Mean local HBA (sc)", 
		xaxt = "n", yaxt = "n", type = "n", main = main)
axis(1, at = 1:13, labels = LETTERS[1:13], lwd = 3)
axis(2, lwd = 3)
for(g in 1:9){
	lines(x = 1:13, y = Chereji.meanLocalHBA[,g], lwd = 3, col = COL[g])
	legend(x = 11.2, y = 1.8 - 0.12*g, legend = paste("sg", g, sep = ""), 
		lwd = 3, col = COL[g], bty = "n")
}
box(lwd = 3)
main <- "Brogaard local HBA (sc)"
plot(x = 1:13, y = numeric(length = 13), ylim = ylim, ylab = "Mean local HBA (sc)", 
		xaxt = "n", yaxt = "n", type = "n", main = main)
axis(1, at = 1:13, labels = LETTERS[1:13], lwd = 3)
axis(2, lwd = 3)
for(g in 1:9){
	lines(x = 1:13, y = Brogaard.meanLocalHBA[,g], lwd = 3, col = COL[g])
	legend(x = 11.2, y = 1.8 - 0.12*g, legend = paste("sg", g, sep = ""), 
		lwd = 3, col = COL[g], bty = "n")
}
box(lwd = 3)
ylim <- c(-0.5801816, 1.6358595)
main <- "Chereji local HBA (sc)"
plot(x = 1:13, y = numeric(length = 13), ylim = ylim, ylab = "Mean local HBA (sc)", 
		xaxt = "n", yaxt = "n", type = "n", main = main)
axis(1, at = 1:13, labels = LETTERS[1:13], lwd = 3)
axis(2, lwd = 3)
for(g in 1:9){
	lines(x = 1:13, y = Chereji.meanLocalHBA[,g], lwd = 3, col = COL[g])
	legend(x = 11.2, y = 1.8 - 0.12*g, legend = paste("sg", g, sep = ""), 
		lwd = 3, col = COL[g], bty = "n")
}
box(lwd = 3)
main <- "Brogaard local HBA (sc)"
plot(x = 1:13, y = numeric(length = 13), ylim = ylim, ylab = "Mean local HBA (sc)", 
		xaxt = "n", yaxt = "n", type = "n", main = main)
axis(1, at = 1:13, labels = LETTERS[1:13], lwd = 3)
axis(2, lwd = 3)
for(g in 1:9){
	lines(x = 1:13, y = Brogaard.meanLocalHBA[,g], lwd = 3, col = COL[g])
	legend(x = 11.2, y = 1.8 - 0.12*g, legend = paste("sg", g, sep = ""), 
		lwd = 3, col = COL[g], bty = "n")
}
box(lwd = 3)


