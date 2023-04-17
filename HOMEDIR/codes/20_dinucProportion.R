#########################################################
## 20_dinucProportion.R
## Dinucleotide proportion analysis 
## (Figure 6 and Supplementary Figure S8): 20210818_nuCpos2_2.txt
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
load(file = "Chereji_w107_add6.RData")
load(file = "Brogaard_Uniq_add6.RData")
library(Biostrings)
setwd(HOMEDIR)
setwd("20210706_parameters")
setwd("RData")
load(file = "merged33_67352_147_nDSS.RData")
load(file = "sc_genome_2011.RData")
genomicFreq <- dinucleotideFrequency(sc.genome.2011, simplify.as = "collapsed", as.prob = TRUE)
Chereji_w107.plus1F <- subset(Chereji_w107, P1 == TRUE & totalNuc >= 2 & strand == "+")
Chereji_w107.plus1R <- subset(Chereji_w107, P1 == TRUE & totalNuc >= 2 & strand == "-")
Chereji_w107.plus1R$seq <- as.character(reverseComplement(DNAStringSet(Chereji_w107.plus1R$seq)))
Chereji_w107.plus1 <- DNAStringSet(rbind(Chereji_w107.plus1F, Chereji_w107.plus1R)$seq)
Chereji_w107.plus2F <- subset(Chereji_w107, P2 == TRUE & totalNuc >= 3 & strand == "+")
Chereji_w107.plus2R <- subset(Chereji_w107, P2 == TRUE & totalNuc >= 3 & strand == "-")
Chereji_w107.plus2R$seq <- as.character(reverseComplement(DNAStringSet(Chereji_w107.plus2R$seq)))
Chereji_w107.plus2 <- DNAStringSet(rbind(Chereji_w107.plus2F, Chereji_w107.plus2R)$seq)
Chereji_w107.plus3F <- subset(Chereji_w107, P3 == TRUE & totalNuc >= 4 & strand == "+")
Chereji_w107.plus3R <- subset(Chereji_w107, P3 == TRUE & totalNuc >= 4 & strand == "-")
Chereji_w107.plus3R$seq <- as.character(reverseComplement(DNAStringSet(Chereji_w107.plus3R$seq)))
Chereji_w107.plus3 <- DNAStringSet(rbind(Chereji_w107.plus3F, Chereji_w107.plus3R)$seq)
Chereji_w107.plus4F <- subset(Chereji_w107, P4 == TRUE & totalNuc >= 5 & strand == "+")
Chereji_w107.plus4R <- subset(Chereji_w107, P4 == TRUE & totalNuc >= 5 & strand == "-")
Chereji_w107.plus4R$seq <- as.character(reverseComplement(DNAStringSet(Chereji_w107.plus4R$seq)))
Chereji_w107.plus4 <- DNAStringSet(rbind(Chereji_w107.plus4F, Chereji_w107.plus4R)$seq)
Brogaard_Uniq.plus1F <- subset(Brogaard_Uniq, P1 == TRUE & totalNuc >= 2 & strand == "+")
Brogaard_Uniq.plus1R <- subset(Brogaard_Uniq, P1 == TRUE & totalNuc >= 2 & strand == "-")
Brogaard_Uniq.plus1R$seq <- as.character(reverseComplement(DNAStringSet(Brogaard_Uniq.plus1R$seq)))
Brogaard_Uniq.plus1 <- DNAStringSet(rbind(Brogaard_Uniq.plus1F, Brogaard_Uniq.plus1R)$seq)
Brogaard_Uniq.plus2F <- subset(Brogaard_Uniq, P2 == TRUE & totalNuc >= 3 & strand == "+")
Brogaard_Uniq.plus2R <- subset(Brogaard_Uniq, P2 == TRUE & totalNuc >= 3 & strand == "-")
Brogaard_Uniq.plus2R$seq <- as.character(reverseComplement(DNAStringSet(Brogaard_Uniq.plus2R$seq)))
Brogaard_Uniq.plus2 <- DNAStringSet(rbind(Brogaard_Uniq.plus2F, Brogaard_Uniq.plus2R)$seq)
Brogaard_Uniq.plus3F <- subset(Brogaard_Uniq, P3 == TRUE & totalNuc >= 4 & strand == "+")
Brogaard_Uniq.plus3R <- subset(Brogaard_Uniq, P3 == TRUE & totalNuc >= 4 & strand == "-")
Brogaard_Uniq.plus3R$seq <- as.character(reverseComplement(DNAStringSet(Brogaard_Uniq.plus3R$seq)))
Brogaard_Uniq.plus3 <- DNAStringSet(rbind(Brogaard_Uniq.plus3F, Brogaard_Uniq.plus3R)$seq)
Brogaard_Uniq.plus4F <- subset(Brogaard_Uniq, P4 == TRUE & totalNuc >= 5 & strand == "+")
Brogaard_Uniq.plus4R <- subset(Brogaard_Uniq, P4 == TRUE & totalNuc >= 5 & strand == "-")
Brogaard_Uniq.plus4R$seq <- as.character(reverseComplement(DNAStringSet(Brogaard_Uniq.plus4R$seq)))
Brogaard_Uniq.plus4 <- DNAStringSet(rbind(Brogaard_Uniq.plus4F, Brogaard_Uniq.plus4R)$seq)
Chereji_w107.plus1.RC <- reverseComplement(Chereji_w107.plus1)
Chereji_w107.plus2.RC <- reverseComplement(Chereji_w107.plus2)
Chereji_w107.plus3.RC <- reverseComplement(Chereji_w107.plus3)
Chereji_w107.plus4.RC <- reverseComplement(Chereji_w107.plus4)
Brogaard_Uniq.plus1.RC <- reverseComplement(Brogaard_Uniq.plus1)
Brogaard_Uniq.plus2.RC <- reverseComplement(Brogaard_Uniq.plus2)
Brogaard_Uniq.plus3.RC <- reverseComplement(Brogaard_Uniq.plus3)
Brogaard_Uniq.plus4.RC <- reverseComplement(Brogaard_Uniq.plus4)
Chereji_w107.no <- DNAStringSet(subset(Chereji_w107, strand == "")$seq)
Brogaard_Uniq.no <- DNAStringSet(subset(Brogaard_Uniq, strand == "")$seq)
Chereji_w107.no.RC <- reverseComplement(Chereji_w107.no)
Brogaard_Uniq.no.RC <- reverseComplement(Brogaard_Uniq.no)
Chereji_w107.no.Combined <- c(Chereji_w107.no, Chereji_w107.no.RC)
Brogaard_Uniq.no.Combined <- c(Brogaard_Uniq.no, Brogaard_Uniq.no.RC)
objNames <- c("merged33_67352.147.nDSS", 
		"Chereji_w107.plus1", "Chereji_w107.plus2", 
		"Chereji_w107.plus3", "Chereji_w107.plus4", "Chereji_w107.no", 
		"Chereji_w107.plus1.RC", "Chereji_w107.plus2.RC", 
		"Chereji_w107.plus3.RC", "Chereji_w107.plus4.RC", 
		"Chereji_w107.no.RC", "Chereji_w107.no.Combined", 
		"Brogaard_Uniq.plus1", "Brogaard_Uniq.plus2", 
		"Brogaard_Uniq.plus3", "Brogaard_Uniq.plus4", "Brogaard_Uniq.no", 
		"Brogaard_Uniq.plus1.RC", "Brogaard_Uniq.plus2.RC", 
		"Brogaard_Uniq.plus3.RC", "Brogaard_Uniq.plus4.RC", 
		"Brogaard_Uniq.no.RC", "Brogaard_Uniq.no.Combined")
WINDOW <- 73
dinuc <- list()
for(i in 1:length(objNames)){
	cat(objNames[i], ", ", sep = "")
	dinuc[[objNames[i]]] <- list()
	dinuc[[objNames[i]]][["dinucFreq"]] <- data.frame(pos = -WINDOW:(WINDOW-1), stringsAsFactors = FALSE)
	for(n in 1:(WINDOW*2+1-1)){
		TARGET <- subseq(get(objNames[i])[width(get(objNames[i])) == WINDOW*2+1], start = n, end = n + 1)
		temp <- dinucleotideFrequency(TARGET, simplify.as = "collapsed", as.prob = TRUE)
	
		for(g in 1:16){
			dinuc[[objNames[i]]][["dinucFreq"]][[names(genomicFreq)[g]]][n] <- temp[g]
		}
	}
	dinuc[[objNames[i]]][["dinucOrderFreq"]] <- data.frame(pos = -WINDOW:(WINDOW-1), stringsAsFactors = FALSE)
	dinuc[[objNames[i]]][["dinucOrderNames"]] <- data.frame(pos = -WINDOW:(WINDOW-1), stringsAsFactors = FALSE)
	for(g in 1:16){
		colName <- paste("Rank", g, sep = "")
		dinuc[[objNames[i]]][["dinucOrderFreq"]][[colName]] <- as.numeric(NA)
		dinuc[[objNames[i]]][["dinucOrderNames"]][[colName]] <- as.character(NA)
	}
	for(n in 1:(WINDOW*2+1-1)){
		temp <- dinuc[[objNames[i]]][["dinucFreq"]][n, 2:17]
		dinucOrderFreq <- temp[order(temp, decreasing = TRUE)]
		dinucOrderNames <- names(dinucOrderFreq)
		dinuc[[objNames[i]]][["dinucOrderFreq"]][n, 2:17] <- dinucOrderFreq
		dinuc[[objNames[i]]][["dinucOrderNames"]][n, 2:17] <- dinucOrderNames
	}
}
setwd(HOMEDIR)
dir.create("20210818_dinuc1")
setwd("20210818_dinuc1")
save(dinuc, file = "dinuc.RData")
setwd(HOMEDIR)
setwd("20210818_dinuc1")
load(file = "dinuc.RData")
setwd(HOMEDIR)
setwd("20210818_dinuc1")
library(RColorBrewer)
COL <- character(length = 16)
for(i in 1:16){
	COL[i] <- rgb(r = 1-i/16, g = 0.1, b = i/16, alpha = 0.8)
}
for(i in 1:length(objNames)){
fileName <- paste(sub(pattern = "\\.", replacement = "_", objNames[i]), "_order_half.pdf", sep = "")
junction.window	<- 3
junctionL <- c((-16-junction.window):(-16+junction.window-1))
junctionR <- c((16-junction.window):(16+junction.window-1))
junctionLR <- c(junctionL, junctionR)
if(i >= 2 & i <= 12){
	cleavage.window <- 10
	cleavageL <- c((-25-cleavage.window):(-25+cleavage.window-1))
	cleavageR <- c((25-cleavage.window):(25+cleavage.window-1))
	cleavageLR <- c(cleavageL, cleavageR)
}
if(i >= 13 & i <= 23){
	cleavage.window <- 10
	cleavageL <- c((-5-cleavage.window):(-5+cleavage.window-1))
	cleavageR <- c((5-cleavage.window):(5+cleavage.window-1))
	cleavageLR <- c(cleavageL, cleavageR)
}
pdf(file = fileName, height = 6, width = 6)
ylim <- c(0, 1)
main <- objNames[i]
plot(x = -73:72, y = numeric(length = 146), type = "n", main = main, ylim = ylim, xlim = c(-73, 0), 
		ylab = "Frequency", 
		xlab = "Start position of dinucleotide relative to dyad (bp)", 
		xaxt = "n", yaxt = "n")
if(i > 1){
	polygon(x = c(min(cleavageL), min(cleavageL), max(cleavageL), max(cleavageL)), 
		y = c(-0.03, -0.01, -0.01, -0.03), col = "gray90", border = NA)
}
polygon(x = c(min(junctionL), min(junctionL), max(junctionL), max(junctionL)), 
	y = c(-0.03, -0.01, -0.01, -0.03), col = "gray80", border = NA)
axis(2, lwd = 2)
axis(1, at = seq(-70, 0, 10), labels = NA, lwd = 2)
axis(1, at = seq(-76, 0, 2), labels = NA, lwd = 1)
if(i == 1)	axis(1, at = seq(-60, 0, 20), labels = c("±60", "±40", "±20", "0"), lwd = 2)
if(i >= 2 & i <= 6)	axis(1, at = seq(-60, 0, 20), labels = c("-60", "-40", "-20", "0"), lwd = 2)
if(i >= 7 & i <= 11)	axis(1, at = seq(-60, 0, 20), labels = c("+60", "+40", "+20", "0"), lwd = 2)
if(i >= 12 & i <= 16)	axis(1, at = seq(-60, 0, 20), labels = c("-60", "-40", "-20", "0"), lwd = 2)
if(i >= 17 & i <= 21)	axis(1, at = seq(-60, 0, 20), labels = c("+60", "+40", "+20", "0"), lwd = 2)
box(lwd = 2)
for(n in 1:(WINDOW*2+1-1)){
	dinucOrderFreq <- dinuc[[objNames[i]]][["dinucOrderFreq"]][n, 2:17]
	y2 <- 1
	for(g in 1:16){
		y1 <- y2
		y2 <- y1 - dinucOrderFreq[g]
		polygon(x = c(n-74-0.5, n-74-0.5, n-74+0.5, n-74+0.5), 
			y = c(y1, y2, y2, y1), border = NA, col = COL[g])
	}
}
abline(h = seq(0, 1, 0.1), lty = 2, lwd = 1)
dev.off()
}
for(i in 1:length(objNames)){
fileName <- paste(sub(pattern = "\\.", replacement = "_", objNames[i]), "_order_half_AATT.pdf", sep = "")
junction.window	<- 3
junctionL <- c((-16-junction.window):(-16+junction.window-1))
junctionR <- c((16-junction.window):(16+junction.window-1))
junctionLR <- c(junctionL, junctionR)
if(i >= 2 & i <= 12){
	cleavage.window <- 10
	cleavageL <- c((-25-cleavage.window):(-25+cleavage.window-1))
	cleavageR <- c((25-cleavage.window):(25+cleavage.window-1))
	cleavageLR <- c(cleavageL, cleavageR)
}
if(i >= 13 & i <= 23){
	cleavage.window <- 10
	cleavageL <- c((-5-cleavage.window):(-5+cleavage.window-1))
	cleavageR <- c((5-cleavage.window):(5+cleavage.window-1))
	cleavageLR <- c(cleavageL, cleavageR)
}
pdf(file = fileName, height = 6, width = 6)
ylim <- c(0, 1)
main <- objNames[i]
plot(x = -73:72, y = numeric(length = 146), type = "n", main = main, ylim = ylim, xlim = c(-73, 0), 
		ylab = "Frequency", 
		xlab = "Start position of dinucleotide relative to dyad (bp)", 
		xaxt = "n", yaxt = "n")
if(i > 1){
	polygon(x = c(min(cleavageL), min(cleavageL), max(cleavageL), max(cleavageL)), 
		y = c(-0.03, -0.01, -0.01, -0.03), col = "gray90", border = NA)
}
polygon(x = c(min(junctionL), min(junctionL), max(junctionL), max(junctionL)), 
	y = c(-0.03, -0.01, -0.01, -0.03), col = "gray80", border = NA)
axis(2, lwd = 2)
axis(1, at = seq(-70, 0, 10), labels = NA, lwd = 2)
axis(1, at = seq(-76, 0, 2), labels = NA, lwd = 1)
if(i == 1)	axis(1, at = seq(-60, 0, 20), labels = c("±60", "±40", "±20", "0"), lwd = 2)
if(i >= 2 & i <= 6)	axis(1, at = seq(-60, 0, 20), labels = c("-60", "-40", "-20", "0"), lwd = 2)
if(i >= 7 & i <= 11)	axis(1, at = seq(-60, 0, 20), labels = c("+60", "+40", "+20", "0"), lwd = 2)
if(i >= 12 & i <= 16)	axis(1, at = seq(-60, 0, 20), labels = c("-60", "-40", "-20", "0"), lwd = 2)
if(i >= 17 & i <= 21)	axis(1, at = seq(-60, 0, 20), labels = c("+60", "+40", "+20", "0"), lwd = 2)
box(lwd = 2)
for(n in 1:(WINDOW*2+1-1)){
	dinucOrderFreq <- dinuc[[objNames[i]]][["dinucOrderFreq"]][n, 2:17]
	dinucOrderNames <- dinuc[[objNames[i]]][["dinucOrderNames"]][n, 2:17]
	y2 <- 1
	for(g in 1:16){
		y1 <- y2
		y2 <- y1 - dinucOrderFreq[g]
		polygon(x = c(n-74-0.5, n-74-0.5, n-74+0.5, n-74+0.5), 
			y = c(y1, y2, y2, y1), border = NA, col = COL[g])
		if(dinucOrderNames[g] == "AA"){
			polygon(x = c(n-74-0.5, n-74-0.5, n-74+0.5, n-74+0.5), 
				y = c(y1, y2, y2, y1), border = "lightgreen", col = "lightgreen", density = 20)
		}
		if(dinucOrderNames[g] == "TT"){
			polygon(x = c(n-74-0.5, n-74-0.5, n-74+0.5, n-74+0.5), 
				y = c(y1, y2, y2, y1), border = "lightpink", col = "lightpink", density = 20)
		}
	}
}
abline(h = seq(0, 1, 0.1), lty = 2, lwd = 1)
dev.off()
}
for(i in 1:length(objNames)){
fileName <- paste(sub(pattern = "\\.", replacement = "_", objNames[i]), "_order_half_dinuc.pdf", sep = "")
junction.window	<- 3
junctionL <- c((-16-junction.window):(-16+junction.window-1))
junctionR <- c((16-junction.window):(16+junction.window-1))
junctionLR <- c(junctionL, junctionR)
if(i >= 2 & i <= 12){
	cleavage.window <- 10
	cleavageL <- c((-25-cleavage.window):(-25+cleavage.window-1))
	cleavageR <- c((25-cleavage.window):(25+cleavage.window-1))
	cleavageLR <- c(cleavageL, cleavageR)
}
if(i >= 13 & i <= 23){
	cleavage.window <- 10
	cleavageL <- c((-5-cleavage.window):(-5+cleavage.window-1))
	cleavageR <- c((5-cleavage.window):(5+cleavage.window-1))
	cleavageLR <- c(cleavageL, cleavageR)
}
pdf(file = fileName, height = 6, width = 6)
ylim <- c(0, 1)
main <- objNames[i]
plot(x = -73:72, y = numeric(length = 146), type = "n", main = main, ylim = ylim, xlim = c(-73, 0), 
		ylab = "Frequency", 
		xlab = "Start position of dinucleotide relative to dyad (bp)", 
		xaxt = "n", yaxt = "n")
if(i > 1){
	polygon(x = c(min(cleavageL), min(cleavageL), max(cleavageL), max(cleavageL)), 
		y = c(-0.03, -0.01, -0.01, -0.03), col = "gray90", border = NA)
}
polygon(x = c(min(junctionL), min(junctionL), max(junctionL), max(junctionL)), 
	y = c(-0.03, -0.01, -0.01, -0.03), col = "gray80", border = NA)
axis(2, lwd = 2)
axis(1, at = seq(-70, 0, 10), labels = NA, lwd = 2)
axis(1, at = seq(-76, 0, 2), labels = NA, lwd = 1)
if(i == 1)	axis(1, at = seq(-60, 0, 20), labels = c("±60", "±40", "±20", "0"), lwd = 2)
if(i >= 2 & i <= 6)	axis(1, at = seq(-60, 0, 20), labels = c("-60", "-40", "-20", "0"), lwd = 2)
if(i >= 7 & i <= 11)	axis(1, at = seq(-60, 0, 20), labels = c("+60", "+40", "+20", "0"), lwd = 2)
if(i >= 12 & i <= 16)	axis(1, at = seq(-60, 0, 20), labels = c("-60", "-40", "-20", "0"), lwd = 2)
if(i >= 17 & i <= 21)	axis(1, at = seq(-60, 0, 20), labels = c("+60", "+40", "+20", "0"), lwd = 2)
box(lwd = 2)
for(n in 1:(WINDOW*2+1-1)){
	dinucFreq <- dinuc[[objNames[i]]][["dinucFreq"]][n, 2:17]
	y2 <- 1
	for(g in 1:16){
		y1 <- y2
		y2 <- y1 - dinucFreq[g]
		polygon(x = c(n-74-0.5, n-74-0.5, n-74+0.5, n-74+0.5), 
			y = c(y1, y2, y2, y1), border = NA, col = COL[g])
	}
}
abline(h = seq(0, 1, 0.1), lty = 2, lwd = 1)
dev.off()
}
for(i in 1:length(objNames)){
fileName <- paste(sub(pattern = "\\.", replacement = "_", objNames[i]), "_order_full.pdf", sep = "")
junction.window	<- 3
junctionL <- c((-16-junction.window):(-16+junction.window-1))
junctionR <- c((16-junction.window):(16+junction.window-1))
junctionLR <- c(junctionL, junctionR)
if(i >= 2 & i <= 12){
	cleavage.window <- 10
	cleavageL <- c((-25-cleavage.window):(-25+cleavage.window-1))
	cleavageR <- c((25-cleavage.window):(25+cleavage.window-1))
	cleavageLR <- c(cleavageL, cleavageR)
}
if(i >= 13 & i <= 23){
	cleavage.window <- 10
	cleavageL <- c((-5-cleavage.window):(-5+cleavage.window-1))
	cleavageR <- c((5-cleavage.window):(5+cleavage.window-1))
	cleavageLR <- c(cleavageL, cleavageR)
}
pdf(file = fileName, height = 6, width = 6)
ylim <- c(0, 1)
main <- objNames[i]
plot(x = -73:72, y = numeric(length = 146), type = "n", main = main, ylim = ylim, xlim = c(-73, 73), 
		ylab = "Frequency", 
		xlab = "Start position of dinucleotide relative to dyad (bp)", 
		xaxt = "n", yaxt = "n")
if(i > 1){
	polygon(x = c(min(cleavageL), min(cleavageL), max(cleavageL), max(cleavageL)), 
		y = c(-0.03, -0.01, -0.01, -0.03), col = "gray90", border = NA)
	polygon(x = c(min(cleavageR), min(cleavageR), max(cleavageR), max(cleavageR)), 
		y = c(-0.03, -0.01, -0.01, -0.03), col = "gray90", border = NA)
}
polygon(x = c(min(junctionL), min(junctionL), max(junctionL), max(junctionL)), 
	y = c(-0.03, -0.01, -0.01, -0.03), col = "gray80", border = NA)
polygon(x = c(min(junctionR), min(junctionR), max(junctionR), max(junctionR)), 
	y = c(-0.03, -0.01, -0.01, -0.03), col = "gray80", border = NA)
axis(2, lwd = 2)
axis(1, at = seq(-70, 70, 10), labels = NA, lwd = 2)
axis(1, at = seq(-74, 74, 2), labels = NA, lwd = 1)
axis(1, at = seq(-60, 60, 20), lwd = 2)
box(lwd = 2)
for(n in 1:(WINDOW*2+1-1)){
	dinucOrderFreq <- dinuc[[objNames[i]]][["dinucOrderFreq"]][n, 2:17]
	y2 <- 1
	for(g in 1:16){
		y1 <- y2
		y2 <- y1 - dinucOrderFreq[g]
		polygon(x = c(n-74-0.5, n-74-0.5, n-74+0.5, n-74+0.5), 
			y = c(y1, y2, y2, y1), border = NA, col = COL[g])
	}
}
abline(h = seq(0, 1, 0.1), lty = 2, lwd = 1)
dev.off()
}
for(i in 1:length(objNames)){
fileName <- paste(sub(pattern = "\\.", replacement = "_", objNames[i]), "_order_full_AATT.pdf", sep = "")
junction.window	<- 3
junctionL <- c((-16-junction.window):(-16+junction.window-1))
junctionR <- c((16-junction.window):(16+junction.window-1))
junctionLR <- c(junctionL, junctionR)
if(i >= 2 & i <= 12){
	cleavage.window <- 10
	cleavageL <- c((-25-cleavage.window):(-25+cleavage.window-1))
	cleavageR <- c((25-cleavage.window):(25+cleavage.window-1))
	cleavageLR <- c(cleavageL, cleavageR)
}
if(i >= 13 & i <= 23){
	cleavage.window <- 10
	cleavageL <- c((-5-cleavage.window):(-5+cleavage.window-1))
	cleavageR <- c((5-cleavage.window):(5+cleavage.window-1))
	cleavageLR <- c(cleavageL, cleavageR)
}
pdf(file = fileName, height = 6, width = 6)
ylim <- c(0, 1)
main <- objNames[i]
plot(x = -73:72, y = numeric(length = 146), type = "n", main = main, ylim = ylim, xlim = c(-73, 73), 
		ylab = "Frequency", 
		xlab = "Start position of dinucleotide relative to dyad (bp)", 
		xaxt = "n", yaxt = "n")
if(i > 1){
	polygon(x = c(min(cleavageL), min(cleavageL), max(cleavageL), max(cleavageL)), 
		y = c(-0.03, -0.01, -0.01, -0.03), col = "gray90", border = NA)
	polygon(x = c(min(cleavageR), min(cleavageR), max(cleavageR), max(cleavageR)), 
		y = c(-0.03, -0.01, -0.01, -0.03), col = "gray90", border = NA)
}
polygon(x = c(min(junctionL), min(junctionL), max(junctionL), max(junctionL)), 
	y = c(-0.03, -0.01, -0.01, -0.03), col = "gray80", border = NA)
polygon(x = c(min(junctionR), min(junctionR), max(junctionR), max(junctionR)), 
	y = c(-0.03, -0.01, -0.01, -0.03), col = "gray80", border = NA)
axis(2, lwd = 2)
axis(1, at = seq(-70, 70, 10), labels = NA, lwd = 2)
axis(1, at = seq(-74, 74, 2), labels = NA, lwd = 1)
axis(1, at = seq(-60, 60, 20), lwd = 2)
box(lwd = 2)
for(n in 1:(WINDOW*2+1-1)){
	dinucOrderFreq <- dinuc[[objNames[i]]][["dinucOrderFreq"]][n, 2:17]
	dinucOrderNames <- dinuc[[objNames[i]]][["dinucOrderNames"]][n, 2:17]
	y2 <- 1
	for(g in 1:16){
		y1 <- y2
		y2 <- y1 - dinucOrderFreq[g]
		polygon(x = c(n-74-0.5, n-74-0.5, n-74+0.5, n-74+0.5), 
			y = c(y1, y2, y2, y1), border = NA, col = COL[g])
		if(dinucOrderNames[g] == "AA"){
			polygon(x = c(n-74-0.5, n-74-0.5, n-74+0.5, n-74+0.5), 
				y = c(y1, y2, y2, y1), border = "lightgreen", col = "lightgreen", density = 20)
		}
		if(dinucOrderNames[g] == "TT"){
			polygon(x = c(n-74-0.5, n-74-0.5, n-74+0.5, n-74+0.5), 
				y = c(y1, y2, y2, y1), border = "lightpink", col = "lightpink", density = 20)
		}
	}
}
abline(h = seq(0, 1, 0.1), lty = 2, lwd = 1)
dev.off()
}
for(i in 1:length(objNames)){
fileName <- paste(sub(pattern = "\\.", replacement = "_", objNames[i]), "_order_full_ATTA.pdf", sep = "")
junction.window	<- 3
junctionL <- c((-16-junction.window):(-16+junction.window-1))
junctionR <- c((16-junction.window):(16+junction.window-1))
junctionLR <- c(junctionL, junctionR)
if(i >= 2 & i <= 12){
	cleavage.window <- 10
	cleavageL <- c((-25-cleavage.window):(-25+cleavage.window-1))
	cleavageR <- c((25-cleavage.window):(25+cleavage.window-1))
	cleavageLR <- c(cleavageL, cleavageR)
}
if(i >= 13 & i <= 23){
	cleavage.window <- 10
	cleavageL <- c((-5-cleavage.window):(-5+cleavage.window-1))
	cleavageR <- c((5-cleavage.window):(5+cleavage.window-1))
	cleavageLR <- c(cleavageL, cleavageR)
}
pdf(file = fileName, height = 6, width = 6)
# par(mfrow = c(4, 2))
ylim <- c(0, 1)
main <- objNames[i]
plot(x = -73:72, y = numeric(length = 146), type = "n", main = main, ylim = ylim, xlim = c(-73, 73), 
		ylab = "Frequency", 
		xlab = "Start position of dinucleotide relative to dyad (bp)", 
		xaxt = "n", yaxt = "n")
if(i > 1){
	polygon(x = c(min(cleavageL), min(cleavageL), max(cleavageL), max(cleavageL)), 
		y = c(-0.03, -0.01, -0.01, -0.03), col = "gray90", border = NA)
	polygon(x = c(min(cleavageR), min(cleavageR), max(cleavageR), max(cleavageR)), 
		y = c(-0.03, -0.01, -0.01, -0.03), col = "gray90", border = NA)
}
polygon(x = c(min(junctionL), min(junctionL), max(junctionL), max(junctionL)), 
	y = c(-0.03, -0.01, -0.01, -0.03), col = "gray80", border = NA)
polygon(x = c(min(junctionR), min(junctionR), max(junctionR), max(junctionR)), 
	y = c(-0.03, -0.01, -0.01, -0.03), col = "gray80", border = NA)
axis(2, lwd = 2)
axis(1, at = seq(-70, 70, 10), labels = NA, lwd = 2)
axis(1, at = seq(-74, 74, 2), labels = NA, lwd = 1)
axis(1, at = seq(-60, 60, 20), lwd = 2)
box(lwd = 2)
for(n in 1:(WINDOW*2+1-1)){
	dinucOrderFreq <- dinuc[[objNames[i]]][["dinucOrderFreq"]][n, 2:17]
	dinucOrderNames <- dinuc[[objNames[i]]][["dinucOrderNames"]][n, 2:17]
	y2 <- 1
	for(g in 1:16){
		y1 <- y2
		y2 <- y1 - dinucOrderFreq[g]
		polygon(x = c(n-74-0.5, n-74-0.5, n-74+0.5, n-74+0.5), 
			y = c(y1, y2, y2, y1), border = NA, col = COL[g])
		if(dinucOrderNames[g] == "AT"){
			polygon(x = c(n-74-0.5, n-74-0.5, n-74+0.5, n-74+0.5), 
				y = c(y1, y2, y2, y1), border = "lightblue", col = "lightblue", density = 20)
		}
		if(dinucOrderNames[g] == "TA"){
			polygon(x = c(n-74-0.5, n-74-0.5, n-74+0.5, n-74+0.5), 
				y = c(y1, y2, y2, y1), border = "white", col = "white", density = 20)
		}
	}
}
abline(h = seq(0, 1, 0.1), lty = 2, lwd = 1)
dev.off()
}


