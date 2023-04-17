#########################################################
## 29_violine.R
## Violine plots (Figure S5): 20220308_nuCpos2_1.txt
library(Biostrings)
library(parallel)
library(vioplot)
library(openxlsx)
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
load(file = "Chereji_w107_add6.RData")
Chereji_w107.noN <- subset(Chereji_w107, strand == "")
Chereji_w107.noN.RC <- Chereji_w107.noN
Chereji_w107.noN.RC$seq <- as.character(reverseComplement(DNAStringSet(Chereji_w107.noN$seq)))
Chereji_w107.noN.sym <- rbind(Chereji_w107.noN, Chereji_w107.noN.RC)
checkNucPattern <- function(pos = 9, SEQ, pattern = "AA", dyadIndex = 74){
	seqStart <- dyadIndex + pos
	seqEnd <- dyadIndex + pos + nchar(pattern) - 1
	if(subseq(SEQ, start = seqStart, end = seqEnd) == pattern) return(TRUE)
	if(subseq(SEQ, start = seqStart, end = seqEnd) != pattern) return(FALSE)
}
setwd(HOMEDIR)
OUTDIR1 <- "20220308_pattern_vs_score"
dir.create(OUTDIR1)
setwd(OUTDIR1)
nucName <- "Chereji_w107.noN.sym"
FWpattern <- "AA"
RCpattern <- as.character(reverseComplement(DNAString(FWpattern)))
OUTDIR2 <- paste(FWpattern, RCpattern, sep = "")
setwd(HOMEDIR)
setwd(OUTDIR1)
if(length(grep(pattern = OUTDIR2, list.files())) == 0) dir.create(OUTDIR2)
setwd(OUTDIR2)
dir.create("violin_plot")
testNuc <- get(nucName)
cat(nucName, " (n=", nrow(testNuc), ")\n", sep = "")
setwd("violin_plot")
positions <- c(-47, 69)
for(p in positions){
	FWindex <- unlist(mcMap(f = checkNucPattern, pos = p, 
			SEQ = testNuc$seq, pattern = FWpattern, dyadIndex = 74, mc.cores = 24), 
			use.names = FALSE)
	RCindex <- unlist(mcMap(f = checkNucPattern, pos = p, 
			SEQ = testNuc$seq, pattern = RCpattern, dyadIndex = 74, mc.cores = 24), 
			use.names = FALSE)
	x <- testNuc[FWindex,]$score
	y <- testNuc[RCindex,]$score
	Welch <- t.test(log10(x), log10(y), method = "Welch")
	if(p < 0)	POS <- paste("m", p, sep = "")
	if(p == 0)	POS <- paste(p, sep = "")
	if(p > 0)	POS <- paste("p", p, sep = "")
	pvalue <- formatC(Welch$p.value, digits = 2, format = "e")
	pvalue2 <- sub(pattern = "\\.", replacement = "_", pvalue)
	fileName <- paste(nucName, "_", POS, "_significant_", pvalue2, ".pdf", sep = "")
	pdf(file = fileName, width = 3, height = 9, pointsize = 15)
	subtitle <- nucName
	maincol <- "black"
	if(Welch$p.value < 0.05) maincol <- "red"
	nameL <- FWpattern
	nameR <- RCpattern
	main <- paste("pos=", p, sep = "")
	ymin <- min(log10(x), log10(y))
	ymax <- max(log10(x), log10(y))
	yrange <- ymax - ymin
	plot(x = c(0.2, 2.8), y = c(NA, NA), 
		ylim = c(ymin - 0.2*yrange, ymax + 0.2*yrange), 
		xaxt = "n", yaxt = "n", xlab = NA, ylab = "Dyad score (log10)", lwd = 2)
	abline(h = median(log10(x)), col = "green", lwd = 2)
	abline(h = median(log10(y)), col = "magenta", lwd = 2)
	vioplot(x = list(log10(x), log10(y)), 
		names = c(nameL, nameR), main = main,
		sub = subtitle, col.main = maincol, add = TRUE)
	points(x = 1, y = median(log10(x)), lwd = 1, bg = "green", cex = 1, pch = 21)
	points(x = 2, y = median(log10(y)), lwd = 1, bg = "magenta", cex = 1, pch = 21)
	box(lwd = 2)
	plineHeight <- ymax + 0.08*yrange
	lines(x = c(1, 2), y = c(plineHeight, plineHeight), lwd = 2)
	text(x = 1.5, y = ymax + 0.13*yrange, labels = paste("p=", pvalue, sep = ""))
	nHeight <- ymin - 0.1*yrange
	text(x = 1, y = nHeight, labels = "n= ", cex = 0.9)
	text(x = 2, y = nHeight, labels = "n= ", cex = 0.9)
	nHeight <- ymin - 0.15*yrange
	text(x = 1, y = nHeight, labels = length(x), cex = 0.9)
	text(x = 2, y = nHeight, labels = length(y), cex = 0.9)
	axis(1, at = c(1, 2), labels = c(FWpattern, RCpattern), lwd = 2)
	axis(2, at = seq(floor(ymin), ceiling(ymax), 0.5), lwd = 2)
	dev.off()
}


