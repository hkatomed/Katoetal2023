#########################################################
## 26_cleavAtSpecificSites.R
## Analysis of cleavage at specific sites 
## (Figure 1E and Supplementary Figure S2): 20211206_nuCpos2_4.txt
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
load(file = "Chereji_w107_add7.RData")
setwd(HOMEDIR)
setwd("20210721_GSE51949")
setwd("RData")
load(file = "Henikoff.RData")
library(RColorBrewer)
library(vioplot)
COL <- brewer.pal(9, "Spectral")
violin.cleavage <- function(NucSet = Chereji_w107, segment = "A", 
				strand = "Watson", clvSite = -2, WINDOW = 0, ylim = c(0, 3.5)){
	temp <- NucSet
	colName <- paste("l", segment, "_scGroup", sep = "")
	Cleavage <- list()
	for(g in 1:9){
		sg <- paste("sg", g, sep = "")
		temp2 <- subset(temp, temp[[colName]] == sg)
		Cleavage[[sg]] <- numeric(length = nrow(temp2))
		for(i in 1:nrow(temp2)){
			chrName <- temp2$chr[i]
			pos <- temp2$pos[i] + clvSite
			Cleavage[[sg]][i] <- sum(Henikoff[[chrName]][[strand]]$score[(pos-WINDOW):(pos+WINDOW)], 
						na.rm = TRUE)
		}
		Cleavage[[sg]] <- Cleavage[[sg]][Cleavage[[sg]] != 0]
		Cleavage[[sg]] <- log10(Cleavage[[sg]])
	}
	ylab <- paste("Cleavage (log10): ", clvSite, "±", WINDOW, sep = "")
	main <- paste("Henikoff cleavage (", colName, ") at ", clvSite, " of ", strand, " strand", sep = "")
	plot(x = integer(length = 9), xlim = c(0.5, 9.5), type = "n", 
			xlab = "Group", ylab = ylab,
			xaxt = "n", yaxt = "n", ylim = ylim, main = main)
	abline(h = median(Cleavage[["sg1"]]), col = "gray", lwd = 3)
	vioplot(x = Cleavage, ylim = ylim, main = main, col = COL, lwd = 3, add = TRUE)
	axis(1, lwd = 3, at = 1:9, labels = names(Cleavage))
	axis(2, lwd = 3, at = seq(0, ylim[2], 1), labels = seq(0, ylim[2], 1))
	axis(2, lwd = 3, at = seq(0, ylim[2], 0.5), labels = NA)
	box(lwd = 3)
}
violin.cleavage(NucSet = Chereji_w107, segment = "A", strand = "Watson", clvSite = -2)
violin.cleavage(NucSet = Chereji_w107, segment = "B", strand = "Watson", clvSite = -2)
violin.cleavage(NucSet = Chereji_w107, segment = "C", strand = "Watson", clvSite = -2)
violin.cleavage(NucSet = Chereji_w107, segment = "D", strand = "Watson", clvSite = -2)
violin.cleavage(NucSet = Chereji_w107, segment = "E", strand = "Watson", clvSite = -2)
violin.cleavage(NucSet = Chereji_w107, segment = "F", strand = "Watson", clvSite = -2)
violin.cleavage(NucSet = Chereji_w107, segment = "G", strand = "Watson", clvSite = -2)
violin.cleavage(NucSet = Chereji_w107, segment = "H", strand = "Watson", clvSite = -2)
violin.cleavage(NucSet = Chereji_w107, segment = "I", strand = "Watson", clvSite = -2)
violin.cleavage(NucSet = Chereji_w107, segment = "J", strand = "Watson", clvSite = -2)
violin.cleavage(NucSet = Chereji_w107, segment = "K", strand = "Watson", clvSite = -2)
violin.cleavage(NucSet = Chereji_w107, segment = "L", strand = "Watson", clvSite = -2)
violin.cleavage(NucSet = Chereji_w107, segment = "M", strand = "Watson", clvSite = -2)
violin.cleavage(NucSet = Chereji_w107, segment = "A", strand = "Crick", clvSite = -5)
violin.cleavage(NucSet = Chereji_w107, segment = "B", strand = "Crick", clvSite = -5)
violin.cleavage(NucSet = Chereji_w107, segment = "C", strand = "Crick", clvSite = -5)
violin.cleavage(NucSet = Chereji_w107, segment = "D", strand = "Crick", clvSite = -5)
violin.cleavage(NucSet = Chereji_w107, segment = "E", strand = "Crick", clvSite = -5)
violin.cleavage(NucSet = Chereji_w107, segment = "F", strand = "Crick", clvSite = -5)
violin.cleavage(NucSet = Chereji_w107, segment = "G", strand = "Crick", clvSite = -5)
violin.cleavage(NucSet = Chereji_w107, segment = "H", strand = "Crick", clvSite = -5)
violin.cleavage(NucSet = Chereji_w107, segment = "I", strand = "Crick", clvSite = -5)
violin.cleavage(NucSet = Chereji_w107, segment = "J", strand = "Crick", clvSite = -5)
violin.cleavage(NucSet = Chereji_w107, segment = "K", strand = "Crick", clvSite = -5)
violin.cleavage(NucSet = Chereji_w107, segment = "L", strand = "Crick", clvSite = -5)
violin.cleavage(NucSet = Chereji_w107, segment = "M", strand = "Crick", clvSite = -5)
violin.cleavage(NucSet = Chereji_w107, segment = "A", strand = "Crick", clvSite = 2)
violin.cleavage(NucSet = Chereji_w107, segment = "B", strand = "Crick", clvSite = 2)
violin.cleavage(NucSet = Chereji_w107, segment = "C", strand = "Crick", clvSite = 2)
violin.cleavage(NucSet = Chereji_w107, segment = "D", strand = "Crick", clvSite = 2)
violin.cleavage(NucSet = Chereji_w107, segment = "E", strand = "Crick", clvSite = 2)
violin.cleavage(NucSet = Chereji_w107, segment = "F", strand = "Crick", clvSite = 2)
violin.cleavage(NucSet = Chereji_w107, segment = "G", strand = "Crick", clvSite = 2)
violin.cleavage(NucSet = Chereji_w107, segment = "H", strand = "Crick", clvSite = 2)
violin.cleavage(NucSet = Chereji_w107, segment = "I", strand = "Crick", clvSite = 2)
violin.cleavage(NucSet = Chereji_w107, segment = "J", strand = "Crick", clvSite = 2)
violin.cleavage(NucSet = Chereji_w107, segment = "K", strand = "Crick", clvSite = 2)
violin.cleavage(NucSet = Chereji_w107, segment = "L", strand = "Crick", clvSite = 2)
violin.cleavage(NucSet = Chereji_w107, segment = "M", strand = "Crick", clvSite = 2)
violin.cleavage(NucSet = Chereji_w107, segment = "A", strand = "Watson", clvSite = 5)
violin.cleavage(NucSet = Chereji_w107, segment = "B", strand = "Watson", clvSite = 5)
violin.cleavage(NucSet = Chereji_w107, segment = "C", strand = "Watson", clvSite = 5)
violin.cleavage(NucSet = Chereji_w107, segment = "D", strand = "Watson", clvSite = 5)
violin.cleavage(NucSet = Chereji_w107, segment = "E", strand = "Watson", clvSite = 5)
violin.cleavage(NucSet = Chereji_w107, segment = "F", strand = "Watson", clvSite = 5)
violin.cleavage(NucSet = Chereji_w107, segment = "G", strand = "Watson", clvSite = 5)
violin.cleavage(NucSet = Chereji_w107, segment = "H", strand = "Watson", clvSite = 5)
violin.cleavage(NucSet = Chereji_w107, segment = "I", strand = "Watson", clvSite = 5)
violin.cleavage(NucSet = Chereji_w107, segment = "J", strand = "Watson", clvSite = 5)
violin.cleavage(NucSet = Chereji_w107, segment = "K", strand = "Watson", clvSite = 5)
violin.cleavage(NucSet = Chereji_w107, segment = "L", strand = "Watson", clvSite = 5)
violin.cleavage(NucSet = Chereji_w107, segment = "M", strand = "Watson", clvSite = 5)


