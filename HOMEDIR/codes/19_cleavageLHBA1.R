#########################################################
## 19_cleavageLHBA1.R
## Analysis of the relationship between H4-S47C–induced cleavage and local HBA
## (Figure 1D and Supplementary Figure S1, FGHDEIJ): 20210721_nuCpos2_10.txt
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
load(file = "Chereji_w107_add5.RData")
load(file = "Brogaard_Uniq_add5.RData")
setwd(HOMEDIR)
setwd("20210721_GSE51949")
setwd("RData")
load(file = "Henikoff.RData")
WINDOW <- 100
temp <- Chereji_w107
colNames <- c("lF_scGroup", "lG_scGroup", "lH_scGroup", "lD_scGroup", "lE_scGroup", "lI_scGroup", "lJ_scGroup")
Cleavage <- list()
for(c in 1:length(colNames)){
	colName <- colNames[c]
	Cleavage[[colName]] <- list()
	Cleavage.Watson <- data.frame(matrix(data = 0, nrow = WINDOW*2 + 1, ncol = 9), stringsAsFactors = FALSE)
	Cleavage.Crick <- data.frame(matrix(data = 0, nrow = WINDOW*2 + 1, ncol = 9), stringsAsFactors = FALSE)
	for(g in 1:9){
		sg <- paste("sg", g, sep = "")
		temp2 <- subset(temp, temp[[colName]] == sg)
		colnames(Cleavage.Watson)[g] <- sg
		colnames(Cleavage.Crick)[g] <- sg
		temp.Watson <- matrix(data = 0, nrow = WINDOW*2 + 1, ncol = nrow(temp2))
		temp.Crick <- matrix(data = 0, nrow = WINDOW*2 + 1, ncol = nrow(temp2))
		for(i in 1:nrow(temp2)){
			chrName <- temp2$chr[i]
			pos <- temp2$pos[i]
			temp.Watson[,i] <- Henikoff[[chrName]][["Watson"]]$score[(pos-WINDOW):(pos+WINDOW)]
			temp.Crick[,i] <- Henikoff[[chrName]][["Crick"]]$score[(pos-WINDOW):(pos+WINDOW)]
		}
		Cleavage.Watson[,g] <- rowMeans(temp.Watson)
		Cleavage.Crick[,g] <- rowMeans(temp.Crick)
	}
	Cleavage[[colName]][["Watson"]] <- Cleavage.Watson
	Cleavage[[colName]][["Crick"]] <- Cleavage.Crick
}
Chereji.Cleavage <- Cleavage
temp <- Brogaard_Uniq
colNames <- c("lF_scGroup", "lG_scGroup", "lH_scGroup", "lD_scGroup", "lE_scGroup", "lI_scGroup", "lJ_scGroup")
Cleavage <- list()
for(c in 1:length(colNames)){
	colName <- colNames[c]
	Cleavage[[colName]] <- list()
	Cleavage.Watson <- data.frame(matrix(data = 0, nrow = WINDOW*2 + 1, ncol = 9), stringsAsFactors = FALSE)
	Cleavage.Crick <- data.frame(matrix(data = 0, nrow = WINDOW*2 + 1, ncol = 9), stringsAsFactors = FALSE)
	for(g in 1:9){
		sg <- paste("sg", g, sep = "")
		temp2 <- subset(temp, temp[[colName]] == sg)
		colnames(Cleavage.Watson)[g] <- sg
		colnames(Cleavage.Crick)[g] <- sg
		temp.Watson <- matrix(data = 0, nrow = WINDOW*2 + 1, ncol = nrow(temp2))
		temp.Crick <- matrix(data = 0, nrow = WINDOW*2 + 1, ncol = nrow(temp2))
		for(i in 1:nrow(temp2)){
			chrName <- temp2$chr[i]
			pos <- temp2$pos[i]
			temp.Watson[,i] <- Henikoff[[chrName]][["Watson"]]$score[(pos-WINDOW):(pos+WINDOW)]
			temp.Crick[,i] <- Henikoff[[chrName]][["Crick"]]$score[(pos-WINDOW):(pos+WINDOW)]
		}
		Cleavage.Watson[,g] <- rowMeans(temp.Watson)
		Cleavage.Crick[,g] <- rowMeans(temp.Crick)
	}
	Cleavage[[colName]][["Watson"]] <- Cleavage.Watson
	Cleavage[[colName]][["Crick"]] <- Cleavage.Crick
}
Brogaard.Cleavage <- Cleavage
library(RColorBrewer)
COL <- brewer.pal(9, "Spectral")
draw.cleavage <- function(colName = "lF_scGroup"){
	ylim <- c(-max(Chereji.Cleavage[[colName]]$Watson, Chereji.Cleavage[[colName]]$Crick, 
			Brogaard.Cleavage[[colName]]$Watson, Brogaard.Cleavage[[colName]]$Crick), 
		max(Chereji.Cleavage[[colName]]$Watson, Chereji.Cleavage[[colName]]$Crick, 
			Brogaard.Cleavage[[colName]]$Watson, Brogaard.Cleavage[[colName]]$Crick))
	main <- paste("Henikoff cleavage (", colName, ")", sep = "")
	plot(x = -100:100, y = numeric(length = 201), ylim = ylim, ylab = "Cleavage (Henikoff)", 
		xaxt = "n", yaxt = "n", type = "n", main = main, xlab = "Distance from dyad (bp)")
	axis(1, lwd = 3)
	axis(2, lwd = 3)
	for(g in 1:9){
		lines(x = -100:100, y = Chereji.Cleavage[[colName]]$Watson[,g], lwd = 3, col = COL[g])
		lines(x = -100:100, y = -Chereji.Cleavage[[colName]]$Crick[,g], lwd = 3, col = COL[g])
		legend(x = 70, y = 88 - 6*g, legend = paste("sg", g, sep = ""), 
			lwd = 3, col = COL[g], bty = "n")
	}
	box(lwd = 3)
}
draw.cleavage(colName = "lF_scGroup")
draw.cleavage(colName = "lG_scGroup")
draw.cleavage(colName = "lH_scGroup")
draw.cleavage(colName = "lD_scGroup")
draw.cleavage(colName = "lE_scGroup")
draw.cleavage(colName = "lI_scGroup")
draw.cleavage(colName = "lJ_scGroup")
draw.cleavage2 <- function(colName = "lF_scGroup"){
	ylim <- c(-80, 80)
	main <- paste("Henikoff cleavage (", colName, ")", sep = "")
	plot(x = -100:100, y = numeric(length = 201), ylim = ylim, ylab = "Cleavage (Henikoff)", 
		xaxt = "n", yaxt = "n", type = "n", main = main, xlab = "Distance from dyad (bp)")
	axis(1, lwd = 3)
	axis(2, lwd = 3)
	abline(h = 0, col = "black", lwd = 3)
	abline(v = 0, col = "lightblue", lwd = 3)
	for(g in 1:9){
		lines(x = -100:100, y = Chereji.Cleavage[[colName]]$Watson[,g], lwd = 3, col = COL[g])
		lines(x = -100:100, y = -Chereji.Cleavage[[colName]]$Crick[,g], lwd = 3, col = COL[g])
		legend(x = 70, y = 88 - 6*g, legend = paste("sg", g, sep = ""), 
			lwd = 3, col = COL[g], bty = "n")
	}
	box(lwd = 3)
}
draw.cleavage2(colName = "lF_scGroup")
draw.cleavage2(colName = "lG_scGroup")
draw.cleavage2(colName = "lH_scGroup")
draw.cleavage2(colName = "lD_scGroup")
draw.cleavage2(colName = "lE_scGroup")
draw.cleavage2(colName = "lI_scGroup")
draw.cleavage2(colName = "lJ_scGroup")


