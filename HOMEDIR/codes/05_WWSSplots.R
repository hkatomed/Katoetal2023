#########################################################
## 05_WWSSplots.R
## WW/SS plots (Figure 3B): 20210708_nuCpos2_2.txt
## Calculation of WW/SS frequency
library(Biostrings)
setwd(HOMEDIR)
setwd("20210706_parameters")
setwd("RData")
load(file = "nature11142_s2_67352_147_nDSS.RData")
load(file = "Chereji_67352_147_nDSS.RData")
load(file = "merged33_67352_147_nDSS.RData")
load(file = "sc_genome_2011.RData")
temp <- dinucleotideFrequency(sc.genome.2011, simplify.as = "collapsed", as.prob = TRUE)
gWW <- sum(temp[c(1,4,13,16)])
gSS <- sum(temp[c(6,7,10,11)])
objNames <- c("nature11142_s2_67352.147.nDSS", "Chereji_67352.147.nDSS", 
		"merged33_67352.147.nDSS")
WINDOW <- 73
WWSS <- list()
for(i in 1:length(objNames)){
	cat(objNames[i], ", ", sep = "")
	WWSS[[objNames[i]]] <- list()
	WWSS[[objNames[i]]][["WW"]] <- numeric(length = WINDOW*2+1-1)
	WWSS[[objNames[i]]][["SS"]] <- numeric(length = WINDOW*2+1-1)
	WWSS[[objNames[i]]][["normWW"]] <- numeric(length = WINDOW*2+1-1)
	WWSS[[objNames[i]]][["normSS"]] <- numeric(length = WINDOW*2+1-1)
	for(n in 1:(WINDOW*2+1-1)){
		TARGET <- subseq(get(objNames[i])[width(get(objNames[i])) == WINDOW*2+1], start = n, end = n + 1)
		temp <- dinucleotideFrequency(TARGET, simplify.as = "collapsed", as.prob = TRUE)
		WW <- sum(temp[c(1,4,13,16)])
		SS <- sum(temp[c(6,7,10,11)])
		normWW <- WW/gWW
		normSS <- SS/gSS
		WWSS[[objNames[i]]][["WW"]][n] <- WW
		WWSS[[objNames[i]]][["SS"]][n] <- SS
		WWSS[[objNames[i]]][["normWW"]][n] <- normWW
		WWSS[[objNames[i]]][["normSS"]][n] <- normSS
	}
}
setwd(HOMEDIR)
setwd("20210706_parameters")
dir.create("20210708_WWSS")
setwd("20210708_WWSS")
save(WWSS, file = "WWSS.RData")
## Draw WW/SS plots
ymax.WW <- as.numeric(NA)
ymax.SS <- as.numeric(NA)
ymax.normWWSS <- as.numeric(NA)
ymin.WW <- as.numeric(NA)
ymin.SS <- as.numeric(NA)
ymin.normWWSS <- as.numeric(NA)
for(i in 1:length(objNames)){
	WW <- WWSS[[objNames[i]]][["WW"]]
	SS <- WWSS[[objNames[i]]][["SS"]]
	normWW <- WWSS[[objNames[i]]][["normWW"]]
	normSS <- WWSS[[objNames[i]]][["normSS"]]
	ymax.WW <- max(ymax.WW, WW, na.rm = TRUE)
	ymin.WW <- min(ymin.WW, WW, na.rm = TRUE)
	ymax.SS <- max(ymax.SS, SS, na.rm = TRUE)
	ymin.SS <- min(ymin.SS, SS, na.rm = TRUE)
	ymax.normWWSS <- max(ymax.normWWSS, normWW, normSS, na.rm = TRUE)
	ymin.normWWSS <- min(ymin.normWWSS, normWW, normSS, na.rm = TRUE)
}
ylim.WW <- c(ymin.WW, ymax.WW)
ylim.SS <- c(ymin.SS, ymax.SS)
ylim.normWWSS <- c(ymin.normWWSS, ymax.normWWSS)
for(i in 1:length(objNames)){
	cat(objNames[i], ", ", sep = "")
	x <- 1:((WINDOW*2+1)-1)
	x <- x - WINDOW -1
	WW <- WWSS[[objNames[i]]][["WW"]]
	SS <- WWSS[[objNames[i]]][["SS"]]
	normWW <- WWSS[[objNames[i]]][["normWW"]]
	normSS <- WWSS[[objNames[i]]][["normSS"]]
	filename <- paste(objNames[i], "_WW_ylimAdjusted.jpeg", sep = "")
	jpeg(file = filename, width = 2000, height = 800, quality = 100, res = 200)
	main <- paste(objNames[i], ", WW", sep = "") 
	plot(x = x, y = numeric(length = length(x)), type = "n", main = main, 
		ylim = ylim.WW, xlim = c(min(x), WINDOW+23), xaxt = "n", yaxt = "n", 
		ylab = "Dinucleotide frequency", xlab = "Center of dinucleotide (bp)")
	lines(x = x + 0.5, y = WW, lwd = 2, col = "red")
	legend(x = max(x), y = ylim.WW[2], legend = "WW", lwd = 2, col = "red", bty = "n")
	axis(side = 1, lwd.ticks = 2)
	axis(side = 2, lwd.ticks = 2)
	box(lwd = 2)
	dev.off()
	filename <- paste(objNames[i], "_SS_ylimAdjusted.jpeg", sep = "")
	jpeg(file = filename, width = 2000, height = 800, quality = 100, res = 200)
	main <- paste(objNames[i], ", SS", sep = "") 
	plot(x = x, y = numeric(length = length(x)), type = "n", main = main, 
		ylim = ylim.SS, xlim = c(min(x), WINDOW+23), xaxt = "n", yaxt = "n", 
		ylab = "Dinucleotide frequency", xlab = "Center of dinucleotide (bp)")
	lines(x = x + 0.5, y = SS, lwd = 2, col = "blue")
	legend(x = max(x), y = ylim.SS[2], legend = "SS", lwd = 2, col = "blue", bty = "n")
	axis(side = 1, lwd.ticks = 2)
	axis(side = 2, lwd.ticks = 2)
	box(lwd = 2)
	dev.off()
	filename <- paste(objNames[i], "_normWWSS_ylimAdjusted.jpeg", sep = "")
	jpeg(file = filename, width = 2000, height = 800, quality = 100, res = 200)
	main <- paste(objNames[i], ", WW/SS (normalized)", sep = "") 
	plot(x = x, y = numeric(length = length(x)), type = "n", main = main, 
		ylim = ylim.normWWSS, xlim = c(min(x), WINDOW+23), xaxt = "n", yaxt = "n", 
		ylab = "Dinucleotide frequency", xlab = "Center of dinucleotide (bp)")
	lines(x = x + 0.5, y = normWW, lwd = 2, col = "red")
	lines(x = x + 0.5, y = normSS, lwd = 2, col = "blue")
	legend(x = max(x), y = ylim.normWWSS[2], legend = "WW", lwd = 2, col = "red", bty = "n")
	legend(x = max(x), y = ylim.normWWSS[2]-0.15, legend = "SS", lwd = 2, col = "blue", bty = "n")
	axis(side = 1, lwd.ticks = 2)
	axis(side = 2, lwd.ticks = 2)
	box(lwd = 2)
	dev.off()
	filename <- paste(objNames[i], "_normWWSS_ylimFixed.jpeg", sep = "")
	jpeg(file = filename, width = 2000, height = 800, quality = 100, res = 200)
	main <- paste(objNames[i], ", WW/SS (normalized)", sep = "") 
	plot(x = x, y = numeric(length = length(x)), type = "n", main = main, 
		ylim = c(0.55, 1.45), xlim = c(min(x), WINDOW+23), xaxt = "n", yaxt = "n", 
		ylab = "Dinucleotide frequency", xlab = "Center of dinucleotide (bp)")
	lines(x = x + 0.5, y = normWW, lwd = 2, col = "red")
	lines(x = x + 0.5, y = normSS, lwd = 2, col = "blue")
	legend(x = max(x), y = ylim.normWWSS[2], legend = "WW", lwd = 2, col = "red", bty = "n")
	legend(x = max(x), y = ylim.normWWSS[2]-0.15, legend = "SS", lwd = 2, col = "blue", bty = "n")
	axis(side = 1, lwd.ticks = 2)
	axis(side = 2, lwd.ticks = 2)
	box(lwd = 2)
	dev.off()
	ylim <- c(min(WW), max(WW))
	filename <- paste(objNames[i], "_WW.jpeg", sep = "")
	jpeg(file = filename, width = 2000, height = 800, quality = 100, res = 200)
	main <- paste(objNames[i], ", WW", sep = "") 
	plot(x = x, y = numeric(length = length(x)), type = "n", main = main, 
		ylim = ylim, xlim = c(min(x), WINDOW+23), xaxt = "n", yaxt = "n", 
		ylab = "Dinucleotide frequency", xlab = "Center of dinucleotide (bp)")
	lines(x = x + 0.5, y = WW, lwd = 2, col = "red")
	legend(x = max(x), y = ylim[2], legend = "WW", lwd = 2, col = "red", bty = "n")
	axis(side = 1, lwd.ticks = 2)
	axis(side = 2, lwd.ticks = 2)
	box(lwd = 2)
	dev.off()
	ylim <- c(min(SS), max(SS))
	filename <- paste(objNames[i], "_SS.jpeg", sep = "")
	jpeg(file = filename, width = 2000, height = 800, quality = 100, res = 200)
	main <- paste(objNames[i], ", SS", sep = "") 
	plot(x = x, y = numeric(length = length(x)), type = "n", main = main, 
		ylim = ylim, xlim = c(min(x), WINDOW+23), xaxt = "n", yaxt = "n", 
		ylab = "Dinucleotide frequency", xlab = "Center of dinucleotide (bp)")
	lines(x = x + 0.5, y = SS, lwd = 2, col = "blue")
	legend(x = max(x), y = ylim[2], legend = "SS", lwd = 2, col = "blue", bty = "n")
	axis(side = 1, lwd.ticks = 2)
	axis(side = 2, lwd.ticks = 2)
	box(lwd = 2)
	dev.off()
	ylim <- c(min(c(normWW, normSS)), max(c(normWW, normSS)))
	filename <- paste(objNames[i], "_normWWSS.jpeg", sep = "")
	jpeg(file = filename, width = 2000, height = 800, quality = 100, res = 200)
	main <- paste(objNames[i], ", WW/SS (normalized)", sep = "") 
	plot(x = x, y = numeric(length = length(x)), type = "n", main = main, 
		ylim = ylim, xlim = c(min(x), WINDOW+23), xaxt = "n", yaxt = "n", 
		ylab = "Dinucleotide frequency", xlab = "Center of dinucleotide (bp)")
	lines(x = x + 0.5, y = normWW, lwd = 2, col = "red")
	lines(x = x + 0.5, y = normSS, lwd = 2, col = "blue")
	legend(x = max(x), y = ylim[2], legend = "WW", lwd = 2, col = "red", bty = "n")
	legend(x = max(x), y = ylim[2]-0.15, legend = "SS", lwd = 2, col = "blue", bty = "n")
	axis(side = 1, lwd.ticks = 2)
	axis(side = 2, lwd.ticks = 2)
	box(lwd = 2)
	dev.off()
}


