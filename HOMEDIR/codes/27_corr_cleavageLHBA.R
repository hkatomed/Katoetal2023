#########################################################
## 27_corr_cleavageLHBA.R
## Peason's correlation test for cleavage and local HBA 
## (Figure 1F and Supplementary Table S1): 20211206_nuCpos2_5.txt
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
clvSites <- c(-5, -2, 2, 5)
strands <- c("Crick", "Watson", "Crick", "Watson")
segments <- LETTERS[1:13]
WINDOW <- 0
temp <- Chereji_w107
df <- data.frame(segment = character(length = 0), strand = character(length = 0), 
			clvSite = integer(length = 0), 
			coefficient = numeric(length = 0), p.value = numeric(length = 0), 
			conf.int.l = numeric(length = 0), conf.int.u = numeric(length = 0))
for(s in 1:length(segments)){
	segment <- segments[s]
	colName <- paste("l", segment, "_sc", sep = "")

	for(c in 1:4){
		clvSite <- clvSites[c]
		strand <- strands[c]
		clvCounts <- numeric(length = nrow(temp))
		for(i in 1:nrow(temp)){
			chrName <- temp$chr[i]
			pos <- temp$pos[i] + clvSite
			clvCounts[i] <- sum(Henikoff[[chrName]][[strand]]$score[(pos-WINDOW):(pos+WINDOW)], 
						na.rm = TRUE)
		}
		x <- temp[[colName]]
		y <- clvCounts
		result <- cor.test(x, y, alternative = "two.sided", method = "pearson", conf.level = 0.95)
		p.value <- result$p.value
		estimate <- result$estimate
		conf.int <- result$conf.int
		df.temp <- data.frame(segment = segment, strand = strand, 
				clvSite = as.integer(clvSite), coefficient = estimate, p.value = p.value, 
				conf.int.l = conf.int[1], conf.int.u = conf.int[2])
		df <- rbind(df, df.temp)
	}
}
library(openxlsx)
setwd(HOMEDIR)
dir.create("20211206_cleavage_correlation")
setwd("20211206_cleavage_correlation")
write.xlsx(df, file = "20211206_cleavage_correlation.xlsx")	# Table S1
save(df, file = "20211206_cleavage_correlation.RData")
clvSites <- c(-5, -2, 2, 5)
ymin <- 1
ymax <- 0
for(c in 1:4){
	df.selected <- subset(df, df$clvSite == clvSites[c])
	ymin <- min(ymin, min(df.selected$conf.int.l))
	ymax <- max(ymax, max(df.selected$conf.int.u))
}
ylim <- c(ymin, ymax)
COLs <- c(	rgb(1, 0.1, 0.2, alpha = 0.8), 
		rgb(1, 0.4, 0.2, alpha = 0.8), 
		rgb(0.2, 0.1, 1, alpha = 0.8), 
		rgb(0.2, 0.4, 1, alpha = 0.8))
plot(x = 1:13, y = numeric(length = 13), ylim = ylim, 
		xaxt = "n", yaxt = "n", type = "n", xlab = "Segment", 
		ylab = "Pearson's correlation coefficient")
axis(1, at = 1:13, labels = LETTERS[1:13], lwd = 3)
axis(2, lwd = 3)
box(lwd = 3)
for(c in 1:4){
	df.selected <- subset(df, df$clvSite == clvSites[c])
	lines(x = 1:13, y = df.selected$coefficient, 
		col = COLs[c], lwd = 3)
	for(s in 1:13){
		arrows(x0 = s, y0 = df.selected$conf.int.l[s], 
			x1 = s, y1 = df.selected$conf.int.u[s], 
			length = 0.05, angle = 90, code = 3, 
			lwd = 3, COLs[c])
	}
	if(clvSites[c] == -5)	Legend <- "-5 Crick"
	if(clvSites[c] == -2)	Legend <- "-2 Watson"
	if(clvSites[c] == 5)	Legend <- "+5 Watson"
	if(clvSites[c] == 2)	Legend <- "+2 Crick"
	legend(x = 9, y = ymax -(c-1)*0.02, legend = Legend, 
			lwd = 3, col = COLs[c], bty = "n")
}


