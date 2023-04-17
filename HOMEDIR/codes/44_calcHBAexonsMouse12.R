#########################################################
## 44_calcHBAexonsMouse12.R
## Calculation of HBA along protein-coding genes (mouse)
## (Figure 4E): 20220309_nuCpos2_1.txt
setwd(HOMEDIR)
setwd("20211210_mouse_GRCm38_p6")
setwd("RData")
load(file = "nucSeq_over500_CentralNuc_HBAlist_CentralNuc.RData")
setwd(HOMEDIR)
setwd("20211210_mouse_GRCm38_p6")
setwd("RData")
load(file = "redNucMatrix_mm10_over500_CentralNuc.RData")
norm.mouse.dyads <- colMeans(redNucMatrix_mm10_over500_CentralNuc)
norm.mouse.dyads <- norm.mouse.dyads/max(norm.mouse.dyads)
setwd(HOMEDIR)
dir.create("20220309_Mouse_HBA")
setwd("20220309_Mouse_HBA")
library(TTR)
scH3H4_33L.HBA <- colMeans(HBAlist_CentralNuc[["chimeric"]][["HBA"]])
scH3H4_33L.LRA1.L <- colMeans(HBAlist_CentralNuc[["chimeric"]][["L"]])
scH3H4_33L.LRA1.R <- colMeans(HBAlist_CentralNuc[["chimeric"]][["R"]])
alpha = 1
alpha.unsmoothed = 0.3
xlim <- c(1, 1001)
pdf.width <- 14
pdf.height <- 10
pdf.pointsize <- 20
ylim <- c(-3.0, 0.5)
fileName <- "mouse_scH3H4_33L_LRA1_dyad_score.pdf"
pdf(file = fileName, width = pdf.width, height = pdf.height, pointsize = pdf.pointsize)
plot(x = 1:1001, y = numeric(length = 1001), type = "n", xaxt = "n", 
	xlab = "Distance from Chereji's +1 dyad (bp)", ylab = "Mean affinity score (log10)", 
	ylim = ylim, xlim = xlim, 
	main = "Model: scH3H4_33L LRA1 and dyad score")
axis(1, at = seq(1, 1001, 100), labels = NA, lwd = 3)
axis(1, at = seq(101, 1001, 200), labels = seq(-400, 400, 200), lwd = 3)
axis(2, lwd = 3)
axis(3, at = seq(1, 1001, 100), labels = NA, lwd = 3)
# lines(x = 1:1001, y = scH3H4_33L.HBA, col = rgb(0, 0.5, 0, alpha = alpha.unsmoothed), lwd = 2)
lines(x = 1:1001, y = scH3H4_33L.LRA1.R, col = rgb(0, 0, 1, alpha = alpha.unsmoothed), lwd = 2)
lines(x = 1:1001, y = scH3H4_33L.LRA1.L, col = rgb(1, 0, 0, alpha = alpha.unsmoothed), lwd = 2)
norm.mouse.dyads <- colMeans(redNucMatrix_mm10_over500_CentralNuc)/
			max(colMeans(redNucMatrix_mm10_over500_CentralNuc))
norm.mouse.dyads <- norm.mouse.dyads * (min(scH3H4_33L.LRA1.R[(501-10):(501+10)]) - ylim[1]) * 8 + ylim[1]
norm.mouse.dyads[501] <- min(scH3H4_33L.LRA1.R[(501-10):(501+10)]) - 0.05
lines(x = 1:1001, y = norm.mouse.dyads, col = gray(0, alpha = 0.9), lwd = 2)
dyads.ymax <- (min(scH3H4_33L.LRA1.R[(501-10):(501+10)]) - ylim[1]) * 8 + ylim[1]
dyads.ymin <- ylim[1]
dyads.scales <- c(dyads.ymin, dyads.ymin+(dyads.ymax-dyads.ymin)*0.05, 
			dyads.ymin+(dyads.ymax-dyads.ymin)*0.10, 
			dyads.ymin+(dyads.ymax-dyads.ymin)*0.15, 
			dyads.ymin+(dyads.ymax-dyads.ymin)*0.20, 
			dyads.ymin+(dyads.ymax-dyads.ymin)*0.25, 
			dyads.ymin+(dyads.ymax-dyads.ymin)*0.30)
axis(4, at = dyads.scales, labels = NA, lwd = 3)
axis(4, at = dyads.scales[c(1,3,5,7)], labels = seq(0, 0.3, 0.1))
temp_SMA.R <- numeric(length = 1001)
temp_SMA.R[1:1001] <- NA
temp_SMA.R[26:(1001-25)] <- SMA(scH3H4_33L.LRA1.R, n = 51)[51:1001]
lines(x = 1:1001, y = temp_SMA.R, type = "l", col = rgb(0, 0, 1, alpha = alpha), lwd = 3)
temp_SMA.L <- numeric(length = 1001)
temp_SMA.L[1:1001] <- NA
temp_SMA.L[26:(1001-25)] <- SMA(scH3H4_33L.LRA1.L, n = 51)[51:1001]
lines(x = 1:1001, y = temp_SMA.L, type = "l", col = rgb(1, 0, 0, alpha = alpha), lwd = 3)
box(lwd = 4)
dev.off()
ylim <- c(-2.5, 0)
alpha = 1
plot(x = 1:1001, y = numeric(length = 1001), type = "n", xaxt = "n", 
	xlab = "Distance from the nucleosomes overlapping exon center (bp)", 
	ylim = ylim, 
	main = "Model: scH3H4_33L HBA vs LRA1")
axis(1, at = seq(1, 1001, 200), labels = seq(-500, 500, 200), lwd = 2)
axis(2, lwd = 2)
temp_SMA <- numeric(length = 1001)
temp_SMA[1:1001] <- NA
temp_SMA[26:(1001-25)] <- SMA(scH3H4_33L.HBA, n = 51)[51:1001]
lines(x = 1:1001, y = temp_SMA, type = "l", col = rgb(0, 0.5, 0, alpha = alpha), lwd = 2)
temp_SMA.R <- numeric(length = 1001)
temp_SMA.R[1:1001] <- NA
temp_SMA.R[26:(1001-25)] <- SMA(scH3H4_33L.LRA1.R, n = 51)[51:1001]
lines(x = 1:1001, y = temp_SMA.R, type = "l", col = rgb(0, 0, 1, alpha = alpha), lwd = 2)
temp_SMA.L <- numeric(length = 1001)
temp_SMA.L[1:1001] <- NA
temp_SMA.L[26:(1001-25)] <- SMA(scH3H4_33L.LRA1.L, n = 51)[51:1001]
lines(x = 1:1001, y = temp_SMA.L, type = "l", col = rgb(1, 0, 0, alpha = alpha), lwd = 2)
legend(x = 750, y = -4, legend = "scH3H4_33L.HBA", col = rgb(0, 0.5, 0, alpha = alpha), lwd = 2, bty = "n")
legend(x = 750, y = -4.4, legend = "scH3H4_33L.LRA1.L", col = rgb(1, 0, 0, alpha = alpha), lwd = 2, bty = "n")
legend(x = 750, y = -4.8, legend = "scH3H4_33L.LRA1.R", col = rgb(0, 0, 1, alpha = alpha), lwd = 2, bty = "n")
box(lwd = 2)
NCP_score <- colMeans(redNucMatrix_mm10_over500_CentralNuc)[1:1001]
temp_SMA.R <- numeric(length = 1001)
temp_SMA.R[1:1001] <- NA
temp_SMA.R2 <- numeric(length = 1001)
temp_SMA.R2[1:1001] <- NA
temp_SMA.R[11:(1001-10)] <- SMA(NCP_score, n = 21)[21:1001]
temp_SMA.R2[21:(1001-20)] <- SMA(temp_SMA.R[11:(1001-10)], n = 21)[21:981]
plot(x = -500:500, 
	y = temp_SMA.R2, 
	type = "l", ylim = c(0.03, 0.2), 
	xaxt = "n", yaxt = "n", ylab = NA, xlab = NA, bty = "n", col = gray(0.1, alpha = 0.8))
axis(4)


