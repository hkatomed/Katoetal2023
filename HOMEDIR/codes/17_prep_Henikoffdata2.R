#########################################################
## 17_prep_Henikoffdata2.R
## Preparation of Henikoff's data: 20210721_nuCpos2_8.txt (2)
## R: Preparation of clevage data for analysis in the R environment
HOMEDIR <- "XXX"
setwd(HOMEDIR)
setwd("20210721_GSE51949/bedGraph")
SAMPLES <- c("SRX372005", "SRX372006", "SRX372007")
for(s in 1:length(SAMPLES)){
	SAMPLE <- SAMPLES[s]
	filename <- paste(SAMPLE, "_Watson.bedGraph", sep = "")
	Watson.bedGraph.each <- read.table(filename, stringsAsFactors = FALSE)
	if(s == 1) Watson.bedGraph <- Watson.bedGraph.each
	if(s != 1) Watson.bedGraph$V3 <- Watson.bedGraph$V3 + Watson.bedGraph.each$V3
	filename <- paste(SAMPLE, "_Crick.bedGraph", sep = "")
	Crick.bedGraph.each <- read.table(filename, stringsAsFactors = FALSE)
	if(s == 1) Crick.bedGraph <- Crick.bedGraph.each
	if(s != 1) Crick.bedGraph$V3 <- Crick.bedGraph$V3 + Crick.bedGraph.each$V3
}
chrNames <- unique(Watson.bedGraph$V1)
Henikoff <- list()
for(c in 1:length(chrNames)){
	chrName <- chrNames[c]
	Henikoff[[chrName]] <- list()
	temp <- subset(Watson.bedGraph, V1 == chrName)[,c(2,3)]
	colnames(temp) <- c("pos", "score")
	temp$score[1:(nrow(temp)-1)] <- temp$score[2:nrow(temp)]
	temp$score[nrow(temp)] <- 0
	Henikoff[[chrName]][["Watson"]] <- temp
	temp <- subset(Crick.bedGraph, V1 == chrName)[,c(2,3)]
	colnames(temp) <- c("pos", "score")
	temp$score[2:nrow(temp)] <- temp$score[1:(nrow(temp)-1)]
	temp$score[nrow(temp)] <- 1
	Henikoff[[chrName]][["Crick"]] <- temp
}
setwd(HOMEDIR)
setwd("20210721_GSE51949")
dir.create("RData")
setwd("RData")
save(Henikoff, file = "Henikoff.RData")


