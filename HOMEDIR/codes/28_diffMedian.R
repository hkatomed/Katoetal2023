#########################################################
## 28_diffMedian.R
## Checking difference of median between groups
## (Table S2): 20211207_nuCpos2_1.txt
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
load(file = "Chereji_w107_add7.RData")
setwd(HOMEDIR)
setwd("20210721_GSE51949")
setwd("RData")
load(file = "Henikoff.RData")
clvSites <- c(-5, -2, 2, 5)
strands <- c("Crick", "Watson", "Crick", "Watson")
segments <- LETTERS[1:13]
WINDOW <- 0
temp <- Chereji_w107
segment <- "A"
strand <- "Watson"
clvSite <- -2
Cleavage <- list()
for(s in 1:length(segments)){
	segment <- segments[s]
	Cleavage[[segment]] <- list()
	colName <- paste("l", segment, "_scGroup", sep = "")
	for(g in 1:9){
		sg <- paste("sg", g, sep = "")
		temp2 <- subset(temp, temp[[colName]] == sg)
		Cleavage[[segment]][[sg]] <- list()
		for(c in 1:4){
			clvSite <- clvSites[c]
			strand <- strands[c]
			if(clvSite < 0) siteName <- paste(strand, "_", clvSite, sep = "")
			if(clvSite > 0) siteName <- paste(strand, "_+", clvSite, sep = "")
			Cleavage[[segment]][[sg]][[siteName]] <- numeric(length = nrow(temp2))
				for(i in 1:nrow(temp2)){
					chrName <- temp2$chr[i]
					pos <- temp2$pos[i] + clvSite
					Cleavage[[segment]][[sg]][[siteName]][i] <- 
						sum(Henikoff[[chrName]][[strand]]$score[(pos-WINDOW):(pos+WINDOW)], 
						na.rm = TRUE)
				}
		}
	}
}
colNames <- names(Cleavage[["A"]][["sg1"]])
for(c in 1:length(colNames)){
	colName <- colNames[c]
	df.each <- data.frame(clvSite = colName, stringsAsFactors = FALSE)
	for(s in 1:length(segments)){
		segment <- segments[s]
		sg1Median <- median(Cleavage[[segment]][["sg1"]][[colName]], na.rm = TRUE)
		sg9Median <- median(Cleavage[[segment]][["sg9"]][[colName]], na.rm = TRUE)
		sg1vs9 <- sg1Median/sg9Median
		df.each[[segment]] <- sg1vs9
	}
	if(c == 1)	df <- df.each
	if(c != 1)	df <- rbind(df, df.each)
}
library(openxlsx)
setwd(HOMEDIR)
dir.create("20211207_median_cleavage_fold")
setwd("20211207_median_cleavage_fold")
write.xlsx(df, file = "20211207_median_cleavage_fold.xlsx")	# Table S2
save(df, file = "20211207_median_cleavage_fold.RData")


