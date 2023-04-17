#########################################################
## 02_prep_H3Q85C.R
## Preparation of the H3-Q85C dataset: 20210706_nuCpos2_1.txt (2)
## R: Preparation of a normalized H3-Q85C dataset.
HOMEDIR <- "XXX"
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
wig <- list()
SAMPLES <- c("GSM2561057", "GSM2561058", "GSM2561059")
for(s in 1:3){
	SAMPLE <- SAMPLES[s]
	filename <- paste(SAMPLE, "_Dyads_H3_CC_rep_", s, ".wig", sep = "")
	wigLines <- readLines(con = filename)
	headers <- grep(pattern = "chr", x = wigLines)

	wig[[SAMPLE]] <- list()
	for(c in 1:16){
		chrName1 <- paste("chr", c, sep = "")
		chrName2 <- paste("chr", as.roman(c), sep = "")
		chrLen <- length(sc.genome.2011[[chrName1]])
		PATTERN <- paste("chrom=", chrName2, " ", sep = "")
		chrHeaderPositions <- grep(pattern = PATTERN, x = wigLines)
		chrHeaders <- wigLines[chrHeaderPositions]
		chr.wig <- data.frame(pos = 1:chrLen, score = as.integer(NA))
		for(h in 1:length(chrHeaderPositions)){
			wigLines.startPos <- chrHeaderPositions[h] + 1
			wigLines.endPos <- chrHeaderPositions[h+1] - 1
			headerNow <- chrHeaders[h]
			headerNext <- chrHeaders[h+1]
			temp <- strsplit(headerNow, split = " ")[[1]][3]
			acctual.startPos <- as.integer(strsplit(temp, split = "start=")[[1]][2])
			temp <- strsplit(headerNext, split = " ")[[1]][3]
			acctual.endPos <- as.integer(strsplit(temp, split = "start=")[[1]][2]) - 1
			if(h == length(chrHeaderPositions)){
				wigLines.endPos <- wigLines.startPos + (chrLen - acctual.startPos)
				acctual.endPos <- chrLen
			}
			scores <- as.numeric(wigLines[wigLines.startPos:wigLines.endPos])
			chr.wig$score[acctual.startPos:acctual.endPos] <- scores
		}
		wig[[SAMPLE]][[chrName1]] <- chr.wig
	}
}
wig[["combined"]] <- list()
for(c in 1:16){
	chrName1 <- paste("chr", c, sep = "")
	chr.wig1 <- wig[[SAMPLES[1]]][[chrName1]]
	chr.wig2 <- wig[[SAMPLES[2]]][[chrName1]]
	chr.wig3 <- wig[[SAMPLES[3]]][[chrName1]]
	chr.wigc <- chr.wig1
	chr.wigc$score <- chr.wig1$score + chr.wig2$score + chr.wig3$score
	wig[["combined"]][[chrName1]] <- chr.wigc
}
wig[["combined2"]] <- list()
for(c in 1:16){
	chrName1 <- paste("chr", c, sep = "")
	chr.wig1 <- wig[[SAMPLES[1]]][[chrName1]]
	chr.wig2 <- wig[[SAMPLES[2]]][[chrName1]]
	chr.wig3 <- wig[[SAMPLES[3]]][[chrName1]]
	chr.wigc <- chr.wig1
	chr.wigc$score <- (chr.wig1$score + chr.wig2$score + chr.wig3$score) /3
	wig[["combined2"]][[chrName1]] <- chr.wigc
}
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
dir.create("RData")
setwd("RData")
save(wig, file = "Chereji_wig.RData")


