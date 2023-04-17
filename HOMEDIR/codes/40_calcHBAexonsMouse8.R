#########################################################
## 40_calcHBAexonsMouse8.R
## Calculation of HBA along protein-coding genes (mouse)
## Preparation of exons: 20211217_nuCpos2_1.txt		matchedExon_over500_CentralNuc.RData
library(Biostrings)
library(parallel)
setwd(HOMEDIR)
setwd("20211210_mouse_GRCm38_p6")
setwd("RData")
load(file = "uniqueNuc_mm10.RData")
load(file = "matchedExon.RData")
matchedExon_over500 <- subset(matchedExon, length >500)
unique(matchedExon_over500$chrNameUCSC)
getNucAtCenter <- function(chrName, exonLeft, exonRight, strand){
	center <- floor((exonLeft + exonRight)/2)
	start <- center - 72
	end <- center + 72
	hit <- subset(uniqueNuc_mm10, chr == chrName & pos >= start & pos <= end)$pos
	if(length(hit) > 1) 	cat("length(hit)>1, i:", i, "\n", sep = "")
	if(length(hit) == 1)	return(hit)
	if(length(hit) == 0)	return(as.integer(NA))
}
matchedExon_over500$centralNuc <- as.integer(NA)
for(i in 1:ceiling(nrow(matchedExon_over500)/1000)){
	from <- (i-1)*1000 + 1
	to <- i*1000
	if(i == ceiling(nrow(matchedExon_over500)/1000))	to <- nrow(matchedExon_over500)
	if(i%%10 == 1) cat("from:", from, " to:", to, "\n", sep = "")
	matchedExon_over500$centralNuc[from:to] <- unlist(mcMap(f = getNucAtCenter, 
				chrName = matchedExon_over500$chrNameUCSC[from:to], 
				exonLeft <- matchedExon_over500$left[from:to], 
				exonRight <- matchedExon_over500$right[from:to], 
				strand <- matchedExon_over500$strand[from:to], 
				mc.cores = 24), use.names = FALSE)
}
for(i in 1:ceiling(nrow(matchedExon_over500)/1000)){
	from <- (i-1)*1000 + 1
	to <- i*1000
	if(i == ceiling(nrow(matchedExon_over500)/1000))	to <- nrow(matchedExon_over500)
	NAsum <- sum(is.na(matchedExon_over500$centralNuc[from:to]))
	nonNAsum <- sum(is.na(matchedExon_over500$centralNuc[from:to])==FALSE)
	cat("# i:", i, " from:", from, " to:", to, " NAsum:", NAsum, 
		" nonNAsum:", nonNAsum, " total:", NAsum+nonNAsum, "\n", sep = "")
}
matchedExon_over500_CentralNuc <- subset(matchedExon_over500, is.na(centralNuc) == FALSE)
setwd(HOMEDIR)
setwd("20211210_mouse_GRCm38_p6")
setwd("RData")
save(matchedExon_over500, file = "matchedExon_over500_withCentralNuc.RData")
save(matchedExon_over500_CentralNuc, file = "matchedExon_over500_CentralNuc.RData")


