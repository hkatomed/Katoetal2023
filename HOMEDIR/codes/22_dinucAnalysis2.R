#########################################################
## 22_dinucAnalysis2.R
## Dinucleotide analysis (RR/YY)
## (Supplementary Figure S7 and Table S3-6): 20211201_nuCpos2_2.txt
library(Biostrings)
library(parallel)
library(vioplot)
library(openxlsx)
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
load(file = "Chereji_w107_add6.RData")
load(file = "Brogaard_Uniq_add6.RData")
Chereji_w107.plus1F <- subset(Chereji_w107, P1 == TRUE & totalNuc >= 2 & strand == "+")
Chereji_w107.plus1R <- subset(Chereji_w107, P1 == TRUE & totalNuc >= 2 & strand == "-")
Chereji_w107.plus1R$seq <- as.character(reverseComplement(DNAStringSet(Chereji_w107.plus1R$seq)))
Chereji_w107.plus1N <- rbind(Chereji_w107.plus1F, Chereji_w107.plus1R)
Chereji_w107.plus2F <- subset(Chereji_w107, P2 == TRUE & totalNuc >= 3 & strand == "+")
Chereji_w107.plus2R <- subset(Chereji_w107, P2 == TRUE & totalNuc >= 3 & strand == "-")
Chereji_w107.plus2R$seq <- as.character(reverseComplement(DNAStringSet(Chereji_w107.plus2R$seq)))
Chereji_w107.plus2N <- rbind(Chereji_w107.plus2F, Chereji_w107.plus2R)
Chereji_w107.plus3F <- subset(Chereji_w107, P3 == TRUE & totalNuc >= 4 & strand == "+")
Chereji_w107.plus3R <- subset(Chereji_w107, P3 == TRUE & totalNuc >= 4 & strand == "-")
Chereji_w107.plus3R$seq <- as.character(reverseComplement(DNAStringSet(Chereji_w107.plus3R$seq)))
Chereji_w107.plus3N <- rbind(Chereji_w107.plus3F, Chereji_w107.plus3R)
Chereji_w107.plus4F <- subset(Chereji_w107, P4 == TRUE & totalNuc >= 5 & strand == "+")
Chereji_w107.plus4R <- subset(Chereji_w107, P4 == TRUE & totalNuc >= 5 & strand == "-")
Chereji_w107.plus4R$seq <- as.character(reverseComplement(DNAStringSet(Chereji_w107.plus4R$seq)))
Chereji_w107.plus4N <- rbind(Chereji_w107.plus4F, Chereji_w107.plus4R)
Brogaard_Uniq.plus1F <- subset(Brogaard_Uniq, P1 == TRUE & totalNuc >= 2 & strand == "+")
Brogaard_Uniq.plus1R <- subset(Brogaard_Uniq, P1 == TRUE & totalNuc >= 2 & strand == "-")
Brogaard_Uniq.plus1R$seq <- as.character(reverseComplement(DNAStringSet(Brogaard_Uniq.plus1R$seq)))
Brogaard_Uniq.plus1N <- rbind(Brogaard_Uniq.plus1F, Brogaard_Uniq.plus1R)
Brogaard_Uniq.plus2F <- subset(Brogaard_Uniq, P2 == TRUE & totalNuc >= 3 & strand == "+")
Brogaard_Uniq.plus2R <- subset(Brogaard_Uniq, P2 == TRUE & totalNuc >= 3 & strand == "-")
Brogaard_Uniq.plus2R$seq <- as.character(reverseComplement(DNAStringSet(Brogaard_Uniq.plus2R$seq)))
Brogaard_Uniq.plus2N <- rbind(Brogaard_Uniq.plus2F, Brogaard_Uniq.plus2R)
Brogaard_Uniq.plus3F <- subset(Brogaard_Uniq, P3 == TRUE & totalNuc >= 4 & strand == "+")
Brogaard_Uniq.plus3R <- subset(Brogaard_Uniq, P3 == TRUE & totalNuc >= 4 & strand == "-")
Brogaard_Uniq.plus3R$seq <- as.character(reverseComplement(DNAStringSet(Brogaard_Uniq.plus3R$seq)))
Brogaard_Uniq.plus3N <- rbind(Brogaard_Uniq.plus3F, Brogaard_Uniq.plus3R)
Brogaard_Uniq.plus4F <- subset(Brogaard_Uniq, P4 == TRUE & totalNuc >= 5 & strand == "+")
Brogaard_Uniq.plus4R <- subset(Brogaard_Uniq, P4 == TRUE & totalNuc >= 5 & strand == "-")
Brogaard_Uniq.plus4R$seq <- as.character(reverseComplement(DNAStringSet(Brogaard_Uniq.plus4R$seq)))
Brogaard_Uniq.plus4N <- rbind(Brogaard_Uniq.plus4F, Brogaard_Uniq.plus4R)
Chereji_w107.noN <- subset(Chereji_w107, strand == "")
Brogaard_Uniq.noN <- subset(Brogaard_Uniq, strand == "")
Chereji_w107.noN.RC <- Chereji_w107.noN
Brogaard_Uniq.noN.RC <- Brogaard_Uniq.noN
Chereji_w107.noN.RC$seq <- as.character(reverseComplement(DNAStringSet(Chereji_w107.noN$seq)))
Brogaard_Uniq.noN.RC$seq <- as.character(reverseComplement(DNAStringSet(Brogaard_Uniq.noN$seq)))
Chereji_w107.noN.sym <- rbind(Chereji_w107.noN, Chereji_w107.noN.RC)
Brogaard_Uniq.noN.sym <- rbind(Brogaard_Uniq.noN, Brogaard_Uniq.noN.RC)
Chereji_w107.RC <- Chereji_w107
Brogaard_Uniq.RC <- Brogaard_Uniq
Chereji_w107.RC$seq <- as.character(reverseComplement(DNAStringSet(Chereji_w107$seq)))
Brogaard_Uniq.RC$seq <- as.character(reverseComplement(DNAStringSet(Brogaard_Uniq$seq)))
Chereji_w107.sym <- rbind(Chereji_w107, Chereji_w107.RC)
Brogaard_Uniq.sym <- rbind(Brogaard_Uniq, Brogaard_Uniq.RC)
nucNames <- c("Chereji_w107.sym", "Brogaard_Uniq.sym", 
		"Chereji_w107.noN.sym", "Brogaard_Uniq.noN.sym", 
		"Chereji_w107.noN", "Chereji_w107.plus1N", "Chereji_w107.plus2N", 
		"Chereji_w107.plus3N", "Chereji_w107.plus4N", 
		"Brogaard_Uniq.noN", "Brogaard_Uniq.plus1N", "Brogaard_Uniq.plus2N", 
		"Brogaard_Uniq.plus3N", "Brogaard_Uniq.plus4N", 
		"Chereji_w107", "Brogaard_Uniq")
checkNucPattern <- function(pos = 9, SEQ, pattern = "AA", dyadIndex = 74){
	seqStart <- dyadIndex + pos
	seqEnd <- dyadIndex + pos + nchar(pattern) - 1
	if(subseq(SEQ, start = seqStart, end = seqEnd) == pattern) return(TRUE)
	if(subseq(SEQ, start = seqStart, end = seqEnd) != pattern) return(FALSE)
}
splitPatterns <- function(pattern){
	if(pattern == "RR"){
		patterns <- c("AA", "AG", "GA", "GG")
	}
	if(pattern == "RRR"){
		patterns <- c("AAA", "AAG", "AGA", "AGG", "GAA", "GAG", "GGA", "GGG")
	}
	return(patterns)
}
setwd(HOMEDIR)
OUTDIR1 <- "20211201_pattern_vs_score2"
dir.create(OUTDIR1)
setwd(OUTDIR1)
PATTERNS <- c("RR", "RRR")
# PATTERNS <- c("RRR")
for(pt in 1:length(PATTERNS)){
	FWpattern <- PATTERNS[pt]
	if(FWpattern == "RR") RCpattern <- "YY"
	if(FWpattern == "RRR") RCpattern <- "YYY"
	FWpatterns <- splitPatterns(pattern = FWpattern)
	RCpatterns <- as.character(reverseComplement(DNAStringSet(FWpatterns)))
	OUTDIR2 <- paste(FWpattern, RCpattern, sep = "")
	setwd(HOMEDIR)
	setwd(OUTDIR1)
	if(length(grep(pattern = OUTDIR2, list.files())) == 0) dir.create(OUTDIR2)
	setwd(OUTDIR2)
	dir.create("xlsx")
	dir.create("RData")
	dir.create("boxplot")
	dir.create("plot")
for(n in 1:length(nucNames)){
	nucName <- nucNames[n]
	testNuc <- get(nucName)
	cat(nucName, " (n=", nrow(testNuc), ")\n", sep = "")
	setwd("boxplot")
	dir.create(nucName)
	setwd("../")
	positions <- c(-73:(73-nchar(FWpattern)+1))
	df1 <- data.frame(nucName = nucName, pos = positions, FW = as.integer(NA), RC = as.integer(NA), 
				FW.proportion = as.numeric(NA), Welch = as.numeric(NA), 
				medianScoreFW = as.integer(NA), medianScoreRC = as.integer(NA), 
				scoreRatioRC = as.numeric(NA))
	for(p in positions){
		for(f in 1:length(FWpatterns)){
			if(f == 1)	tempF <- unlist(mcMap(f = checkNucPattern, pos = p, 
				SEQ = testNuc$seq, pattern = FWpatterns[f], dyadIndex = 74, mc.cores = 4), 
				use.names = FALSE)
			if(f != 1)	tempF <- tempF | unlist(mcMap(f = checkNucPattern, pos = p, 
				SEQ = testNuc$seq, pattern = FWpatterns[f], dyadIndex = 74, mc.cores = 4), 
				use.names = FALSE)
		}
		for(r in 1:length(RCpatterns)){
			if(r == 1)	tempR <- unlist(mcMap(f = checkNucPattern, pos = p, 
				SEQ = testNuc$seq, pattern = RCpatterns[r], dyadIndex = 74, mc.cores = 4), 
				use.names = FALSE)
			if(r != 1)	tempR <- tempR | unlist(mcMap(f = checkNucPattern, pos = p, 
				SEQ = testNuc$seq, pattern = RCpatterns[r], dyadIndex = 74, mc.cores = 4), 
				use.names = FALSE)
		}
		FWindex <- tempF
		RCindex <- tempR
		i <- p + 74
		df1$FW[i] <- nrow(testNuc[FWindex,])
		df1$RC[i] <- nrow(testNuc[RCindex,])
		df1$FW.proportion[i] <- df1$FW[i]/(df1$FW[i] + df1$RC[i])
		x <- testNuc[FWindex,]$score
		y <- testNuc[RCindex,]$score
		df1$medianScoreFW[i] <- median(x, na.rm = TRUE)
		df1$medianScoreRC[i] <- median(y, na.rm = TRUE)
		df1$scoreRatioRC[i] <- df1$medianScoreFW[i]/df1$medianScoreRC[i]
		if(length(x) > 1 & length(y) > 1){
			df1$Welch[i] <- t.test(log10(x), log10(y), method = "Welch")$p.value
			if(df1$Welch[i] < 0.05){
				setwd("boxplot")
				setwd(nucName)
				if(p < 0)	POS <- paste("m", p, sep = "")
				if(p == 0)	POS <- paste(p, sep = "")
				if(p > 0)	POS <- paste("p", p, sep = "")
				pvalue <- formatC(df1$Welch[i], format = "e")
				pvalue2 <- sub(pattern = "\\.", replacement = "_", pvalue)
				fileName <- paste(nucName, "_", POS, "_significant_", pvalue2, ".pdf", sep = "")
				pdf(fileName)
				subtitle <- paste("p=", pvalue, sep = "")
				maincol <- "black"
				if(df1$Welch[i] < 0.05) maincol <- "red"
				nameL <- paste(FWpattern, " (n=", length(x), ")", sep = "")
				nameR <- paste(RCpattern, " (n=", length(y), ")", sep = "")
				main <- paste(nucName, ", ", POS, sep = "")
				vioplot(x = list(log10(x), log10(y)), ylab = "score(log10)",  
					names = c(nameL, nameR), main = main,
					sub = subtitle, col.main = maincol)
				abline(h = mean(log10(x)), col = "red")
				abline(h = mean(log10(y)), col = "blue")
				dev.off()
				setwd("../../")
			}
		}
	}
	setwd("xlsx")
	fileName <- paste(nucName, ".xlsx", sep = "")
	write.xlsx(df1, file = fileName, overwrite = TRUE)
	setwd("../RData")
	fileName <- paste(nucName, ".RData", sep = "")
	save(df1, file = fileName)
	setwd("../")
	setwd("plot")
	fileName <- paste(nucName, "_1.pdf", sep = "")
	pdf(fileName, width = 6, height = 8)
	par(mfrow = c(3,1))
	x <- df1$pos
	if(nchar(FWpattern) == 2){
		x <- x + 0.5
		xlab <- "Center of dinucleotide relative to dyad (bp)"
	}
	if(nchar(FWpattern) == 3){
		x <- x + 1
		xlab <- "Center of trinucleotide relative to dyad (bp)"
	}
	y1 <- df1$FW
	y2 <- df1$RC
	ylab <- "Occurrence"
	ymax <- max(c(y1, y2), na.rm = TRUE)
	ymin <- min(c(y1, y2), na.rm = TRUE)
	ymax <- ymax + 0.2*ymax
	ymin <- ymin - 0.2*ymin
	main <- paste(nucName, " (n=", nrow(testNuc), ")", sep = "")
	if(length(grep(pattern = "sym", nucName)) == 1){
		main <- paste(nucName, " (n=", nrow(testNuc)/2, ")", sep = "")
	}
	plot(x, y = numeric(length = length(y1)), ylim = c(ymin, ymax), 
		type = "n", xaxt = "n", yaxt = "n", xlab = xlab, ylab = ylab, 
		main = main, cex.lab = 1.5)
	legend(x = 54, y = ymax +0.05*(ymax-ymin), legend = FWpattern, lwd = 2, 
			col = rgb(0.2, 0.8, 0.2, alpha = 0.8), bty = "n")
	legend(x = 54, y = ymax -0.05*(ymax-ymin), legend = RCpattern, lwd = 2, 
			col = rgb(0.8, 0.0, 0.5, alpha = 0.8), bty = "n")
	abline(v = 0, lwd = 2, col = "gray")
	lines(x, y1, lwd = 2, col = rgb(0.2, 0.8, 0.2, alpha = 0.8))
	lines(x, y2, lwd = 2, col = rgb(0.8, 0.0, 0.5, alpha = 0.8))
	axis(1, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	y <- df1$FW.proportion
	ylab <- paste(FWpattern, " proportion", sep = "")
	ymax <- max(y, na.rm = TRUE)
	ymin <- min(y, na.rm = TRUE)
	ymax <- ymax + 0.2*ymax
	ymin <- ymin - 0.2*ymin
	plot(x, y, ylim = c(ymin, ymax), 
		type = "n", xaxt = "n", yaxt = "n", xlab = xlab, ylab = ylab, 
		main = paste(nucName, " (n=", nrow(testNuc), ")", sep = ""), cex.lab = 1.5)
	abline(v = 0, lwd = 2, col = "gray")
	abline(h = seq(0, 100, 1), lwd = 2, col = "gray")
	abline(h = 0.5, lwd = 2, col = "orange")	
	lines(x, y, lwd = 2)
	axis(1, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	y <- -log10(df1$Welch)
	ylab <- "-log10 (p-value)"
	plot(x, y, type = "n", xaxt = "n", yaxt = "n", xlab = xlab, ylab = ylab, 
		main = paste(nucName, " (n=", nrow(testNuc), ")", sep = ""), cex.lab = 1.5)
	abline(v = 0, lwd = 2, col = "gray")
	abline(h = -log10(0.05), lwd = 2, col = "orange")
	abline(h = -log10(0.01), lwd = 2, col = "red")
	lines(x, y, lwd = 2)
	axis(1, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	dev.off()
	setwd("../")
	setwd("plot")
	fileName <- paste(nucName, "_2.pdf", sep = "")
	pdf(fileName, width = 6, height = 8)
	par(mfrow = c(3,1))
	y1 <- df1$medianScoreFW
	y2 <- df1$medianScoreRC
	ylab <- "median score"
	ymax <- max(c(y1, y2), na.rm = TRUE)
	ymin <- min(c(y1, y2), na.rm = TRUE)
	ymax <- ymax + 0.3*(ymax-ymin)
	ymin <- ymin - 0.3*(ymax-ymin)
	plot(x, y = numeric(length = length(y1)), ylim = c(ymin, ymax), 
		type = "n", xaxt = "n", yaxt = "n", xlab = xlab, ylab = ylab, 
		main = paste(nucName, " (n=", nrow(testNuc), ")", sep = ""), cex.lab = 1.5)
	legend(x = 54, y = ymax + 0.05*(ymax-ymin), legend = FWpattern, lwd = 2, 
			col = rgb(0.2, 0.8, 0.2, alpha = 0.8), bty = "n")
	legend(x = 54, y = ymax - 0.05*(ymax-ymin), legend = RCpattern, lwd = 2, 
			col = rgb(0.8, 0.0, 0.5, alpha = 0.8), bty = "n")
	abline(v = 0, lwd = 2, col = "gray")
	lines(x, y1, lwd = 2, col = rgb(0.2, 0.8, 0.2, alpha = 0.8))
	lines(x, y2, lwd = 2, col = rgb(0.8, 0.0, 0.5, alpha = 0.8))
	axis(1, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	y <- df1$scoreRatioRC
	ylab <- paste(FWpattern, "/", RCpattern, " ratio (median)", sep = "")
	ymax <- max(y, na.rm = TRUE)
	ymin <- min(y, na.rm = TRUE)
	ymax <- ymax + 0.2*(ymax-ymin)
	ymin <- ymin - 0.2*(ymax-ymin)
	plot(x, y, ylim = c(ymin, ymax), 
		type = "n", xaxt = "n", yaxt = "n", xlab = xlab, ylab = ylab, 
		main = paste(nucName, " (n=", nrow(testNuc), ")", sep = ""), cex.lab = 1.5)
	abline(v = 0, lwd = 2, col = "gray")
	abline(h = seq(0, 100, 1), lwd = 2, col = "gray")
	abline(h = 1, lwd = 2, col = "orange")	
	lines(x, y, lwd = 2)
	axis(1, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	plot(x, y, ylim = c(ymin, ymax), 
		type = "n", xaxt = "n", yaxt = "n", xlab = xlab, ylab = ylab, 
		main = paste(nucName, " (n=", nrow(testNuc), ")", sep = ""), cex.lab = 1.5)
	abline(v = 0, lwd = 2, col = "gray")
	abline(h = seq(0, 100, 1), lwd = 2, col = "gray")
	abline(h = 1, lwd = 2, col = "orange")	
	lines(x, y, lwd = 2)
	axis(1, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	y05 <- y
	y05[df1$Welch > 0.05] <- NA
	y01 <- y
	y01[df1$Welch > 0.01] <- NA
	points(x, y05, col = "white", pch = 20, cex = 0.8)	
	points(x, y05, col = "black", cex = 0.8)		
	points(x, y01, col = "black", pch = 20, cex = 0.8)	
	legend(x = 54, y = ymax + 0.05*(ymax-ymin), legend = "p<0.01", pch = 19, cex = 0.8, 
			col = "black", bty = "n")
	legend(x = 54, y = ymax - 0.05*(ymax-ymin), legend = "p<0.05", pch = 1, cex = 0.8, 
			col = "black", bty = "n")
	dev.off()
	setwd("../")
	setwd("plot")
	fileName <- paste(nucName, "_3.pdf", sep = "")
	pdf(fileName, width = 6, height = 10)
	par(mfrow = c(4,1))
	y1 <- df1$FW
	y2 <- df1$RC
	ylab <- "Occurrence"
	ymax <- max(c(y1, y2), na.rm = TRUE)
	ymin <- min(c(y1, y2), na.rm = TRUE)
	ymax <- ymax + 0.2*ymax
	ymin <- ymin - 0.2*ymin
	main <- paste(nucName, " (n=", nrow(testNuc), ")", sep = "")
	if(length(grep(pattern = "sym", nucName)) == 1){
		main <- paste(nucName, " (n=", nrow(testNuc)/2, ")", sep = "")
	}
	plot(x, y = numeric(length = length(y1)), ylim = c(ymin, ymax), 
		type = "n", xaxt = "n", yaxt = "n", xlab = xlab, ylab = ylab, 
		main = main, cex.lab = 1.5)
	legend(x = 54, y = ymax +0.05*(ymax-ymin), legend = FWpattern, lwd = 2, 
			col = rgb(0.2, 0.8, 0.2, alpha = 0.8), bty = "n")
	legend(x = 54, y = ymax -0.05*(ymax-ymin), legend = RCpattern, lwd = 2, 
			col = rgb(0.8, 0.0, 0.5, alpha = 0.8), bty = "n")
	abline(v = 0, lwd = 2, col = "gray")
	lines(x, y1, lwd = 2, col = rgb(0.2, 0.8, 0.2, alpha = 0.8))
	lines(x, y2, lwd = 2, col = rgb(0.8, 0.0, 0.5, alpha = 0.8))
	axis(1, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	y105 <- y1
	y105[df1$Welch > 0.05] <- NA
	y105[df1$scoreRatioRC < 1] <- NA	
	y101 <- y1
	y101[df1$Welch > 0.01] <- NA
	y101[df1$scoreRatioRC < 1] <- NA	
	points(x, y105, col = "white", pch = 20, cex = 0.8)	
	points(x, y105, col = "black", cex = 0.8)		
	points(x, y101, col = "black", pch = 20, cex = 0.8)	
	y205 <- y2
	y205[df1$Welch > 0.05] <- NA
	y205[df1$scoreRatioRC > 1] <- NA	
	y201 <- y2
	y201[df1$Welch > 0.01] <- NA
	y201[df1$scoreRatioRC > 1] <- NA	
	points(x, y205, col = "white", pch = 20, cex = 0.8)	
	points(x, y205, col = "black", cex = 0.8)		
	points(x, y201, col = "black", pch = 20, cex = 0.8)	
	legend(x = 54, y = ymin + 0.25*(ymax-ymin), legend = "p<0.01", pch = 19, cex = 0.8, 
			col = "black", bty = "n")
	legend(x = 54, y = ymin + 0.15*(ymax-ymin), legend = "p<0.05", pch = 1, cex = 0.8, 
			col = "black", bty = "n")
	y <- df1$FW.proportion
	ylab <- paste(FWpattern, " proportion", sep = "")
	ymax <- max(y, na.rm = TRUE)
	ymin <- min(y, na.rm = TRUE)
	ymax <- ymax + 0.2*ymax
	ymin <- ymin - 0.2*ymin
	plot(x, y, ylim = c(ymin, ymax), 
		type = "n", xaxt = "n", yaxt = "n", xlab = xlab, ylab = ylab, 
		main = paste(nucName, " (n=", nrow(testNuc), ")", sep = ""), cex.lab = 1.5)
	abline(v = 0, lwd = 2, col = "gray")
	abline(h = seq(0, 100, 1), lwd = 2, col = "gray")
	abline(h = 0.5, lwd = 2, col = "orange")	
	lines(x, y, lwd = 2)
	axis(1, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	y05 <- y
	y05[df1$Welch > 0.05] <- NA
	y01 <- y
	y01[df1$Welch > 0.01] <- NA
	points(x, y05, col = "white", pch = 20, cex = 0.8)	
	points(x, y05, col = "black", cex = 0.8)		
	points(x, y01, col = "black", pch = 20, cex = 0.8)	
	legend(x = 54, y = ymin + 0.25*(ymax-ymin), legend = "p<0.01", pch = 19, cex = 0.8, 
			col = "black", bty = "n")
	legend(x = 54, y = ymin + 0.15*(ymax-ymin), legend = "p<0.05", pch = 1, cex = 0.8, 
			col = "black", bty = "n")
	y1 <- df1$medianScoreFW
	y2 <- df1$medianScoreRC
	ylab <- "median score"
	ymax <- max(c(y1, y2), na.rm = TRUE)
	ymin <- min(c(y1, y2), na.rm = TRUE)
	ymax <- ymax + 0.3*(ymax-ymin)
	ymin <- ymin - 0.3*(ymax-ymin)
	plot(x, y = numeric(length = length(y1)), ylim = c(ymin, ymax), 
		type = "n", xaxt = "n", yaxt = "n", xlab = xlab, ylab = ylab, 
		main = paste(nucName, " (n=", nrow(testNuc), ")", sep = ""), cex.lab = 1.5)
	legend(x = 54, y = ymax + 0.05*(ymax-ymin), legend = FWpattern, lwd = 2, 
			col = rgb(0.2, 0.8, 0.2, alpha = 0.8), bty = "n")
	legend(x = 54, y = ymax - 0.05*(ymax-ymin), legend = RCpattern, lwd = 2, 
			col = rgb(0.8, 0.0, 0.5, alpha = 0.8), bty = "n")
	abline(v = 0, lwd = 2, col = "gray")
	lines(x, y1, lwd = 2, col = rgb(0.2, 0.8, 0.2, alpha = 0.8))
	lines(x, y2, lwd = 2, col = rgb(0.8, 0.0, 0.5, alpha = 0.8))
	axis(1, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	y105 <- y1
	y105[df1$Welch > 0.05] <- NA
	y105[df1$scoreRatioRC < 1] <- NA	
	y101 <- y1
	y101[df1$Welch > 0.01] <- NA
	y101[df1$scoreRatioRC < 1] <- NA	
	points(x, y105, col = "white", pch = 20, cex = 0.8)	
	points(x, y105, col = "black", cex = 0.8)		
	points(x, y101, col = "black", pch = 20, cex = 0.8)	
	y205 <- y2
	y205[df1$Welch > 0.05] <- NA
	y205[df1$scoreRatioRC > 1] <- NA	
	y201 <- y2
	y201[df1$Welch > 0.01] <- NA
	y201[df1$scoreRatioRC > 1] <- NA	
	points(x, y205, col = "white", pch = 20, cex = 0.8)	
	points(x, y205, col = "black", cex = 0.8)		
	points(x, y201, col = "black", pch = 20, cex = 0.8)	
	legend(x = 54, y = ymin + 0.25*(ymax-ymin), legend = "p<0.01", pch = 19, cex = 0.8, 
			col = "black", bty = "n")
	legend(x = 54, y = ymin + 0.15*(ymax-ymin), legend = "p<0.05", pch = 1, cex = 0.8, 
			col = "black", bty = "n")
	y <- df1$scoreRatioRC
	ylab <- paste(FWpattern, "/", RCpattern, " ratio (median)", sep = "")
	ymax <- max(y, na.rm = TRUE)
	ymin <- min(y, na.rm = TRUE)
	ymax <- ymax + 0.2*(ymax-ymin)
	ymin <- ymin - 0.2*(ymax-ymin)
	plot(x, y, ylim = c(ymin, ymax), 
		type = "n", xaxt = "n", yaxt = "n", xlab = xlab, ylab = ylab, 
		main = paste(nucName, " (n=", nrow(testNuc), ")", sep = ""), cex.lab = 1.5)
	abline(v = 0, lwd = 2, col = "gray")
	abline(h = seq(0, 100, 1), lwd = 2, col = "gray")
	abline(h = 1, lwd = 2, col = "orange")		
	lines(x, y, lwd = 2)
	axis(1, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	y05 <- y
	y05[df1$Welch > 0.05] <- NA
	y01 <- y
	y01[df1$Welch > 0.01] <- NA
	points(x, y05, col = "white", pch = 20, cex = 0.8)	
	points(x, y05, col = "black", cex = 0.8)		
	points(x, y01, col = "black", pch = 20, cex = 0.8)	
	legend(x = 54, y = ymin + 0.25*(ymax-ymin), legend = "p<0.01", pch = 19, cex = 0.8, 
			col = "black", bty = "n")
	legend(x = 54, y = ymin + 0.15*(ymax-ymin), legend = "p<0.05", pch = 1, cex = 0.8, 
			col = "black", bty = "n")
	dev.off()
	setwd("../")
	setwd("plot")
	fileName <- paste(nucName, "_4.pdf", sep = "")
	pdf(fileName, width = 6, height = 10)
	par(mfrow = c(4,1))
	x <- df1$pos
	xlab <- "Position of dinucleotide relative to dyad (bp)"
	y1 <- df1$FW
	y2 <- df1$RC
	ylab <- "Occurrence"
	ymax <- max(c(y1, y2), na.rm = TRUE)
	ymin <- min(c(y1, y2), na.rm = TRUE)
	ymax <- ymax + 0.2*ymax
	ymin <- ymin - 0.2*ymin
	main <- paste(nucName, " (n=", nrow(testNuc), ")", sep = "")
	if(length(grep(pattern = "sym", nucName)) == 1){
		main <- paste(nucName, " (n=", nrow(testNuc)/2, ")", sep = "")
	}
	plot(x, y = numeric(length = length(y1)), ylim = c(ymin, ymax), 
		type = "n", xaxt = "n", yaxt = "n", xlab = xlab, ylab = ylab, 
		main = main, cex.lab = 1.5)
	legend(x = 54, y = ymax +0.05*(ymax-ymin), legend = FWpattern, lwd = 2, 
			col = rgb(0.2, 0.8, 0.2, alpha = 0.8), bty = "n")
	legend(x = 54, y = ymax -0.05*(ymax-ymin), legend = RCpattern, lwd = 2, 
			col = rgb(0.8, 0.0, 0.5, alpha = 0.8), bty = "n")
	abline(v = 0, lwd = 2, col = "gray")
	lines(x, y1, lwd = 2, col = rgb(0.2, 0.8, 0.2, alpha = 0.8))
	lines(x, y2, lwd = 2, col = rgb(0.8, 0.0, 0.5, alpha = 0.8))
	axis(1, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	y <- df1$FW.proportion
	ylab <- paste(FWpattern, " proportion", sep = "")
	ymax <- max(y, na.rm = TRUE)
	ymin <- min(y, na.rm = TRUE)
	ymax <- ymax + 0.2*ymax
	ymin <- ymin - 0.2*ymin
	plot(x, y, ylim = c(ymin, ymax), 
		type = "n", xaxt = "n", yaxt = "n", xlab = xlab, ylab = ylab, 
		main = paste(nucName, " (n=", nrow(testNuc), ")", sep = ""), cex.lab = 1.5)
	abline(v = 0, lwd = 2, col = "gray")
	abline(h = seq(0, 100, 1), lwd = 2, col = "gray")
	abline(h = 0.5, lwd = 2, col = "orange")	
	lines(x, y, lwd = 2)
	axis(1, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	y1 <- df1$medianScoreFW
	y2 <- df1$medianScoreRC
	ylab <- "median score"
	ymax <- max(c(y1, y2), na.rm = TRUE)
	ymin <- min(c(y1, y2), na.rm = TRUE)
	ymax <- ymax + 0.3*(ymax-ymin)
	ymin <- ymin - 0.3*(ymax-ymin)
	plot(x, y = numeric(length = length(y1)), ylim = c(ymin, ymax), 
		type = "n", xaxt = "n", yaxt = "n", xlab = xlab, ylab = ylab, 
		main = paste(nucName, " (n=", nrow(testNuc), ")", sep = ""), cex.lab = 1.5)
	legend(x = 54, y = ymax + 0.05*(ymax-ymin), legend = FWpattern, lwd = 2, 
			col = rgb(0.2, 0.8, 0.2, alpha = 0.8), bty = "n")
	legend(x = 54, y = ymax - 0.05*(ymax-ymin), legend = RCpattern, lwd = 2, 
			col = rgb(0.8, 0.0, 0.5, alpha = 0.8), bty = "n")
	abline(v = 0, lwd = 2, col = "gray")
	lines(x, y1, lwd = 2, col = rgb(0.2, 0.8, 0.2, alpha = 0.8))
	lines(x, y2, lwd = 2, col = rgb(0.8, 0.0, 0.5, alpha = 0.8))
	axis(1, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	y <- df1$scoreRatioRC
	ylab <- paste(FWpattern, "/", RCpattern, " ratio (median)", sep = "")
	ymax <- max(y, na.rm = TRUE)
	ymin <- min(y, na.rm = TRUE)
	ymax <- ymax + 0.2*(ymax-ymin)
	ymin <- ymin - 0.2*(ymax-ymin)
	plot(x, y, ylim = c(ymin, ymax), 
		type = "n", xaxt = "n", yaxt = "n", xlab = xlab, ylab = ylab, 
		main = paste(nucName, " (n=", nrow(testNuc), ")", sep = ""), cex.lab = 1.5)
	abline(v = 0, lwd = 2, col = "gray")
	abline(h = seq(0, 100, 1), lwd = 2, col = "gray")
	abline(h = 1, lwd = 2, col = "orange")	
	lines(x, y, lwd = 2)
	axis(1, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	dev.off()
	setwd("../")
	setwd("plot")
	fileName <- paste(nucName, "_3_noLegend.pdf", sep = "")
	pdf(fileName, width = 6, height = 10)
	par(mfrow = c(4,1))
	y1 <- df1$FW
	y2 <- df1$RC
	ylab <- "Occurrence"
	ymax <- max(c(y1, y2), na.rm = TRUE)
	ymin <- min(c(y1, y2), na.rm = TRUE)
	ymax <- ymax + 0.2*ymax
	ymin <- ymin - 0.2*ymin
	main <- paste(nucName, " (n=", nrow(testNuc), ")", sep = "")
	if(length(grep(pattern = "sym", nucName)) == 1){
		main <- paste(nucName, " (n=", nrow(testNuc)/2, ")", sep = "")
	}
	plot(x, y = numeric(length = length(y1)), ylim = c(ymin, ymax), 
		type = "n", xaxt = "n", yaxt = "n", xlab = xlab, ylab = ylab, 
		main = main, cex.lab = 1.5)
	abline(v = 0, lwd = 2, col = "gray")
	lines(x, y1, lwd = 2, col = rgb(0.2, 0.8, 0.2, alpha = 0.8))
	lines(x, y2, lwd = 2, col = rgb(0.8, 0.0, 0.5, alpha = 0.8))
	axis(1, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	y105 <- y1
	y105[df1$Welch > 0.05] <- NA
	y105[df1$scoreRatioRC < 1] <- NA	
	y101 <- y1
	y101[df1$Welch > 0.01] <- NA
	y101[df1$scoreRatioRC < 1] <- NA	
	points(x, y105, col = "white", pch = 20, cex = 0.8)		
	points(x, y105, col = "black", cex = 0.8)			
	points(x, y101, col = "black", pch = 20, cex = 0.8)		
	y205 <- y2
	y205[df1$Welch > 0.05] <- NA
	y205[df1$scoreRatioRC > 1] <- NA	
	y201 <- y2
	y201[df1$Welch > 0.01] <- NA
	y201[df1$scoreRatioRC > 1] <- NA	
	points(x, y205, col = "white", pch = 20, cex = 0.8)		
	points(x, y205, col = "black", cex = 0.8)			
	points(x, y201, col = "black", pch = 20, cex = 0.8)
	y <- df1$FW.proportion
	ylab <- paste(FWpattern, " proportion", sep = "")
	ymax <- max(y, na.rm = TRUE)
	ymin <- min(y, na.rm = TRUE)
	ymax <- ymax + 0.2*ymax
	ymin <- ymin - 0.2*ymin
	plot(x, y, ylim = c(ymin, ymax), 
		type = "n", xaxt = "n", yaxt = "n", xlab = xlab, ylab = ylab, 
		main = paste(nucName, " (n=", nrow(testNuc), ")", sep = ""), cex.lab = 1.5)
	abline(v = 0, lwd = 2, col = "gray")
	abline(h = seq(0, 100, 1), lwd = 2, col = "gray")
	abline(h = 0.5, lwd = 2, col = "orange")	
	lines(x, y, lwd = 2)
	axis(1, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	y05 <- y
	y05[df1$Welch > 0.05] <- NA
	y01 <- y
	y01[df1$Welch > 0.01] <- NA
	points(x, y05, col = "white", pch = 20, cex = 0.8)	
	points(x, y05, col = "black", cex = 0.8)		
	points(x, y01, col = "black", pch = 20, cex = 0.8)
	y1 <- df1$medianScoreFW
	y2 <- df1$medianScoreRC
	ylab <- "median score"
	ymax <- max(c(y1, y2), na.rm = TRUE)
	ymin <- min(c(y1, y2), na.rm = TRUE)
	ymax <- ymax + 0.3*(ymax-ymin)
	ymin <- ymin - 0.3*(ymax-ymin)
	plot(x, y = numeric(length = length(y1)), ylim = c(ymin, ymax), 
		type = "n", xaxt = "n", yaxt = "n", xlab = xlab, ylab = ylab, 
		main = paste(nucName, " (n=", nrow(testNuc), ")", sep = ""), cex.lab = 1.5)
	abline(v = 0, lwd = 2, col = "gray")
	lines(x, y1, lwd = 2, col = rgb(0.2, 0.8, 0.2, alpha = 0.8))
	lines(x, y2, lwd = 2, col = rgb(0.8, 0.0, 0.5, alpha = 0.8))
	axis(1, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	y105 <- y1
	y105[df1$Welch > 0.05] <- NA
	y105[df1$scoreRatioRC < 1] <- NA	
	y101 <- y1
	y101[df1$Welch > 0.01] <- NA
	y101[df1$scoreRatioRC < 1] <- NA	
	points(x, y105, col = "white", pch = 20, cex = 0.8)	
	points(x, y105, col = "black", cex = 0.8)		
	points(x, y101, col = "black", pch = 20, cex = 0.8)	
	y205 <- y2
	y205[df1$Welch > 0.05] <- NA
	y205[df1$scoreRatioRC > 1] <- NA	
	y201 <- y2
	y201[df1$Welch > 0.01] <- NA
	y201[df1$scoreRatioRC > 1] <- NA	
	points(x, y205, col = "white", pch = 20, cex = 0.8)	
	points(x, y205, col = "black", cex = 0.8)		
	points(x, y201, col = "black", pch = 20, cex = 0.8)	
	y <- df1$scoreRatioRC
	ylab <- paste(FWpattern, "/", RCpattern, " ratio (median)", sep = "")
	ymax <- max(y, na.rm = TRUE)
	ymin <- min(y, na.rm = TRUE)
	ymax <- ymax + 0.2*(ymax-ymin)
	ymin <- ymin - 0.2*(ymax-ymin)
	plot(x, y, ylim = c(ymin, ymax), 
		type = "n", xaxt = "n", yaxt = "n", xlab = xlab, ylab = ylab, 
		main = paste(nucName, " (n=", nrow(testNuc), ")", sep = ""), cex.lab = 1.5)
	abline(v = 0, lwd = 2, col = "gray")
	abline(h = seq(0, 100, 1), lwd = 2, col = "gray")
	abline(h = 1, lwd = 2, col = "orange")	
	lines(x, y, lwd = 2)
	axis(1, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	y05 <- y
	y05[df1$Welch > 0.05] <- NA
	y01 <- y
	y01[df1$Welch > 0.01] <- NA
	points(x, y05, col = "white", pch = 20, cex = 0.8)	
	points(x, y05, col = "black", cex = 0.8)		
	points(x, y01, col = "black", pch = 20, cex = 0.8)
	dev.off()
	setwd("../")
	rm(df1)
	rm(testNuc)
	gc(reset = TRUE)
} # n-loop
} # pt-loop


