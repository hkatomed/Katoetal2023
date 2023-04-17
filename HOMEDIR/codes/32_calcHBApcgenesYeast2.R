#########################################################
## 32_calcHBApcgenesYeast2.R
## Calculation of HBA along protein-coding genes (yeast)
## (Figure 4B): 20220308_nuCpos2_2.txt
library(parallel)
library(Biostrings)
setwd(HOMEDIR)
setwd("20210906_Chereji_HBA")
setwd("RData")
load(file = "matrix1646_Chereji.RData")	
load(file = "matrix1646.LRA1.L_Chereji.RData")
load(file = "matrix1646.LRA1.R_Chereji.RData")
setwd(HOMEDIR)
dir.create("20220308_Chereji_HBA")
setwd("20220308_Chereji_HBA")
library(TTR)
scH3H4_33L.HBA <- rowMeans(matrix1646, na.rm = TRUE)
scH3H4_33L.LRA1.L <- rowMeans(matrix1646.LRA1.L, na.rm = TRUE)
scH3H4_33L.LRA1.R <- rowMeans(matrix1646.LRA1.R, na.rm = TRUE)
alpha = 1
alpha.unsmoothed = 0.3
xlim <- c(100, 901)
ylim <- c(-5.5, 2)
fileName <- "scH3H4_33L_HBA.pdf"
pdf(file = fileName, width = 8, height = 8)
plot(x = 1:1501, y = numeric(length = 1501), type = "n", xaxt = "n", 
	xlab = "Distance from Chereji's +1 dyad (bp)", ylab = "Mean affinity score (log10)", 
	ylim = ylim, xlim = xlim, 
	main = "Model: scH3H4_33L HBA vs LRA1")
axis(1, at = seq(1, 1501, 100), labels = seq(-500, 1000, 100), lwd = 2)
axis(2, lwd = 2)
lines(x = 1:1501, y = scH3H4_33L.HBA, col = rgb(0, 0.5, 0, alpha = alpha.unsmoothed), lwd = 2)
temp_SMA <- numeric(length = 1501)
temp_SMA[1:1501] <- NA
temp_SMA[26:(1501-25)] <- SMA(scH3H4_33L.HBA, n = 51)[51:1501]
lines(x = 1:1501, y = temp_SMA, type = "l", col = rgb(0, 0.5, 0, alpha = alpha), lwd = 2)
lines(x = c(501, 501), y = c(1.8, 3), col = "gray", lwd = 2)
lines(x = c(501, 501), y = c(-8, -3.3), col = "gray", lwd = 2)
box(lwd = 2)
dev.off()
ylim <- c(-3.5, 1.2)
fileName <- "scH3H4_33L_LRA1.pdf"
pdf(file = fileName, width = 8, height = 8)
plot(x = 1:1501, y = numeric(length = 1501), type = "n", xaxt = "n", 
	xlab = "Distance from Chereji's +1 dyad (bp)", ylab = "Mean affinity score (log10)", 
	ylim = ylim, xlim = xlim, 
	main = "Model: scH3H4_33L HBA vs LRA1")
axis(1, at = seq(1, 1501, 100), labels = seq(-500, 1000, 100), lwd = 2)
axis(2, lwd = 2)
# lines(x = 1:1501, y = scH3H4_33L.HBA, col = rgb(0, 0.5, 0, alpha = alpha.unsmoothed), lwd = 2)
lines(x = 1:1501, y = scH3H4_33L.LRA1.R, col = rgb(0, 0, 1, alpha = alpha.unsmoothed), lwd = 2)
lines(x = 1:1501, y = scH3H4_33L.LRA1.L, col = rgb(1, 0, 0, alpha = alpha.unsmoothed), lwd = 2)
temp_SMA.R <- numeric(length = 1501)
temp_SMA.R[1:1501] <- NA
temp_SMA.R[26:(1501-25)] <- SMA(scH3H4_33L.LRA1.R, n = 51)[51:1501]
lines(x = 1:1501, y = temp_SMA.R, type = "l", col = rgb(0, 0, 1, alpha = alpha), lwd = 2)
temp_SMA.L <- numeric(length = 1501)
temp_SMA.L[1:1501] <- NA
temp_SMA.L[26:(1501-25)] <- SMA(scH3H4_33L.LRA1.L, n = 51)[51:1501]
lines(x = 1:1501, y = temp_SMA.L, type = "l", col = rgb(1, 0, 0, alpha = alpha), lwd = 2)
lines(x = c(501, 501), y = c(1.1, 2), col = "gray", lwd = 2)
lines(x = c(501, 501), y = c(-5, -1.8), col = "gray", lwd = 2)
box(lwd = 2)
dev.off()
ylim <- c(-5.5, 2)
fileName <- "scH3H4_33L_LRA1_and_HBA.pdf"
pdf(file = fileName, width = 8, height = 8)
plot(x = 1:1501, y = numeric(length = 1501), type = "n", xaxt = "n", 
	xlab = "Distance from Chereji's +1 dyad (bp)", ylab = "Mean affinity score (log10)", 
	ylim = ylim, xlim = xlim, 
	main = "Model: scH3H4_33L HBA vs LRA1")
axis(1, at = seq(1, 1501, 100), labels = seq(-500, 1000, 100), lwd = 2)
axis(2, lwd = 2)
lines(x = 1:1501, y = scH3H4_33L.HBA, col = rgb(0, 0.5, 0, alpha = alpha.unsmoothed), lwd = 2)
lines(x = 1:1501, y = scH3H4_33L.LRA1.R, col = rgb(0, 0, 1, alpha = alpha.unsmoothed), lwd = 2)
lines(x = 1:1501, y = scH3H4_33L.LRA1.L, col = rgb(1, 0, 0, alpha = alpha.unsmoothed), lwd = 2)
temp_SMA <- numeric(length = 1501)
temp_SMA[1:1501] <- NA
temp_SMA[26:(1501-25)] <- SMA(scH3H4_33L.HBA, n = 51)[51:1501]
lines(x = 1:1501, y = temp_SMA, type = "l", col = rgb(0, 0.5, 0, alpha = alpha), lwd = 2)
temp_SMA.R <- numeric(length = 1501)
temp_SMA.R[1:1501] <- NA
temp_SMA.R[26:(1501-25)] <- SMA(scH3H4_33L.LRA1.R, n = 51)[51:1501]
lines(x = 1:1501, y = temp_SMA.R, type = "l", col = rgb(0, 0, 1, alpha = alpha), lwd = 2)
temp_SMA.L <- numeric(length = 1501)
temp_SMA.L[1:1501] <- NA
temp_SMA.L[26:(1501-25)] <- SMA(scH3H4_33L.LRA1.L, n = 51)[51:1501]
lines(x = 1:1501, y = temp_SMA.L, type = "l", col = rgb(1, 0, 0, alpha = alpha), lwd = 2)
box(lwd = 2)
dev.off()
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("RData")
load(file = "Chereji_wig.RData")
library(Biostrings)
setwd(HOMEDIR)
setwd("20210706_parameters")
setwd("RData")
load(file = "sc_genome_2011.RData")
setwd(HOMEDIR)
setwd("20210708_Chereji_HBA")
setwd("RData")
load(file = "Chereji_over1500.RData")
Chereji.dyads <- matrix(nrow = 1501, ncol = nrow(Chereji.over1500))
for(i in 1:nrow(Chereji.over1500)){
	chrName <- Chereji.over1500$Chr[i]
	strand <- Chereji.over1500$Strand[i]
	geneName <- Chereji.over1500$ORF[i]
	plus1 <- Chereji.over1500$plus1[i]
	if(strand == "1"){
		Start <- plus1 - 500
		End <- plus1 + 1000
	}
	if(strand == "-1"){
		Start <- plus1 + 500
		End <- plus1 - 1000
	}
	Chereji.dyads[,i] <- wig[["combined"]][[chrName]]$score[Start:End]
}
norm.Chereji.dyads <- rowMeans(Chereji.dyads)
norm.Chereji.dyads <- norm.Chereji.dyads/max(norm.Chereji.dyads)
setwd(HOMEDIR)
dir.create("20220308_Chereji_HBA")
setwd("20220308_Chereji_HBA")
alpha = 1
alpha.unsmoothed = 0.3
xlim <- c(1, 1501)
pdf.width <- 14
pdf.height <- 10
pdf.pointsize <- 20
ylim <- c(-5.5, 2.2)
fileName <- "scH3H4_33L_HBA_dyad_score.pdf"
pdf(file = fileName, width = pdf.width, height = pdf.height, pointsize = pdf.pointsize)
plot(x = 1:1501, y = numeric(length = 1501), type = "n", xaxt = "n", 
	xlab = "Distance from Chereji's +1 dyad (bp)", ylab = "Mean affinity score (log10)", 
	ylim = ylim, xlim = xlim, 
	main = "Model: scH3H4_33L HBA  and dyad score")
axis(1, at = seq(1, 1501, 100), labels = NA, lwd = 3)
axis(1, at = seq(101, 1501, 200), labels = seq(-400, 1000, 200), lwd = 3)
axis(2, lwd = 3)
axis(3, at = seq(1, 1501, 100), labels = NA, lwd = 3)
lines(x = 1:1501, y = scH3H4_33L.HBA, col = rgb(0, 0.5, 0, alpha = alpha.unsmoothed), lwd = 2)
norm.Chereji.dyads <- rowMeans(Chereji.dyads)/max(rowMeans(Chereji.dyads))
norm.Chereji.dyads <- norm.Chereji.dyads * (min(scH3H4_33L.HBA[(501-10):(501+10)]) - ylim[1]) * 0.9 + ylim[1]
lines(x = 1:1501, y = norm.Chereji.dyads, col = gray(0, alpha = 0.9), lwd = 2)
dyads.ymax <- (min(scH3H4_33L.HBA[(501-10):(501+10)]) - ylim[1]) * 0.9 + ylim[1]
dyads.ymin <- ylim[1]
dyads.scales <- c(dyads.ymin, dyads.ymin+(dyads.ymax-dyads.ymin)*0.25, 
			dyads.ymin+(dyads.ymax-dyads.ymin)*0.5, 
			dyads.ymin+(dyads.ymax-dyads.ymin)*0.75, dyads.ymax)
axis(4, at = dyads.scales, labels = NA, lwd = 3)
axis(4, at = dyads.scales[c(1,3,5)], labels = seq(0, 1, 0.5))
temp_SMA <- numeric(length = 1501)
temp_SMA[1:1501] <- NA
temp_SMA[26:(1501-25)] <- SMA(scH3H4_33L.HBA, n = 51)[51:1501]
lines(x = 1:1501, y = temp_SMA, type = "l", col = rgb(0, 0.5, 0, alpha = alpha), lwd = 3)
lines(x = c(673, 673), y = c(-4.7, -1.8), col = "gray", lwd = 2)
lines(x = c(835, 835), y = c(-4.7, -1.8), col = "gray", lwd = 2)
lines(x = c(1003, 1003), y = c(-4.7, -1.8), col = "gray", lwd = 2)
lines(x = c(1167, 1167), y = c(-4.7, -1.8), col = "gray", lwd = 2)
lines(x = c(1328, 1328), y = c(-4.7, -1.8), col = "gray", lwd = 2)
box(lwd = 4)
dev.off()
ylim <- c(-3.5, 1.3)
fileName <- "scH3H4_33L_LRA1_dyad_score.pdf"
pdf(file = fileName, width = pdf.width, height = pdf.height, pointsize = pdf.pointsize)
plot(x = 1:1501, y = numeric(length = 1501), type = "n", xaxt = "n", 
	xlab = "Distance from Chereji's +1 dyad (bp)", ylab = "Mean affinity score (log10)", 
	ylim = ylim, xlim = xlim, 
	main = "Model: scH3H4_33L LRA1 and dyad score")
axis(1, at = seq(1, 1501, 100), labels = NA, lwd = 3)
axis(1, at = seq(101, 1501, 200), labels = seq(-400, 1000, 200), lwd = 3)
axis(2, lwd = 3)
axis(3, at = seq(1, 1501, 100), labels = NA, lwd = 3)
# lines(x = 1:1501, y = scH3H4_33L.HBA, col = rgb(0, 0.5, 0, alpha = alpha.unsmoothed), lwd = 2)
lines(x = 1:1501, y = scH3H4_33L.LRA1.R, col = rgb(0, 0, 1, alpha = alpha.unsmoothed), lwd = 2)
lines(x = 1:1501, y = scH3H4_33L.LRA1.L, col = rgb(1, 0, 0, alpha = alpha.unsmoothed), lwd = 2)
norm.Chereji.dyads <- rowMeans(Chereji.dyads)/max(rowMeans(Chereji.dyads))
norm.Chereji.dyads <- norm.Chereji.dyads * (min(scH3H4_33L.LRA1.R[(501-10):(501+10)]) - ylim[1]) * 0.9 + ylim[1]
lines(x = 1:1501, y = norm.Chereji.dyads, col = gray(0, alpha = 0.9), lwd = 2)
dyads.ymax <- (min(scH3H4_33L.LRA1.R[(501-10):(501+10)]) - ylim[1]) * 0.9 + ylim[1]
dyads.ymin <- ylim[1]
dyads.scales <- c(dyads.ymin, dyads.ymin+(dyads.ymax-dyads.ymin)*0.25, 
			dyads.ymin+(dyads.ymax-dyads.ymin)*0.5, 
			dyads.ymin+(dyads.ymax-dyads.ymin)*0.75, dyads.ymax)
axis(4, at = dyads.scales, labels = NA, lwd = 3)
axis(4, at = dyads.scales[c(1,3,5)], labels = seq(0, 1, 0.5))
temp_SMA.R <- numeric(length = 1501)
temp_SMA.R[1:1501] <- NA
temp_SMA.R[26:(1501-25)] <- SMA(scH3H4_33L.LRA1.R, n = 51)[51:1501]
lines(x = 1:1501, y = temp_SMA.R, type = "l", col = rgb(0, 0, 1, alpha = alpha), lwd = 3)
temp_SMA.L <- numeric(length = 1501)
temp_SMA.L[1:1501] <- NA
temp_SMA.L[26:(1501-25)] <- SMA(scH3H4_33L.LRA1.L, n = 51)[51:1501]
lines(x = 1:1501, y = temp_SMA.L, type = "l", col = rgb(1, 0, 0, alpha = alpha), lwd = 3)
lines(x = c(673, 673), y = c(-3, -1), col = "gray", lwd = 2)
lines(x = c(835, 835), y = c(-3, -1), col = "gray", lwd = 2)
lines(x = c(1003, 1003), y = c(-3, -1), col = "gray", lwd = 2)
lines(x = c(1167, 1167), y = c(-3, -1), col = "gray", lwd = 2)
lines(x = c(1328, 1328), y = c(-3, -1), col = "gray", lwd = 2)
box(lwd = 4)
dev.off()


