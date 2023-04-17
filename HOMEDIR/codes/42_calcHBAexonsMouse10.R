#########################################################
## 42_calcHBAexonsMouse10.R
## Calculation of HBA along protein-coding genes (mouse)
## Preparation of nucleosomes: 20211217_nuCpos2_3.txt	redNucMatrix_mm10_over500_CentralNuc.RData
library(Biostrings)
library(parallel)
setwd(HOMEDIR)
setwd("20211210_mouse_GRCm38_p6")
setwd("RData")
load(file = "matchedExon_over500_CentralNuc.RData")
chrNames <- unique(matchedExon_over500_CentralNuc$chrNameUCSC)
setwd(HOMEDIR)
setwd("GSE82127_RAW")
inDir <- getwd()
redNucList <- list()
for(c in 1:length(chrNames)){
	chrName <- chrNames[c]
	cat(chrName, "\n")
	selectedExons <- subset(matchedExon_over500_CentralNuc, chrNameUCSC == chrNames[c])
	setwd(inDir)
	redundantFile <- paste(chrName, "_Chemical_NCPscore_mm10.txt", sep = "")	
	redundantNuc <- read.table(file = redundantFile, header = FALSE)

	redNucMatrix <- matrix(nrow = nrow(selectedExons), ncol = 1001)
	rownames(redNucMatrix) <- selectedExons$UniqueName
	for(i in 1:nrow(selectedExons)){
		if(i%%10 == 0) cat(i, ",", sep = "")
		centralNuc <- selectedExons$centralNuc[i]
		strand <- selectedExons$strand[i]
		
		left <- centralNuc - 500
		right <- centralNuc + 500
		
		selectedNuc <- subset(redundantNuc, V2 >= left & V2 <= right)
		
		if(strand == "+"){
			selectedNuc$refPos <- selectedNuc$V2 -left + 1
		}
		if(strand == "-"){
			selectedNuc$refPos <- right - selectedNuc$V2 + 1
		}
		nucScores <- integer(length = 1001)
		nucScores[selectedNuc$refPos] <- selectedNuc$V3
		redNucMatrix[i,] <- nucScores
	}
	redNucList[[chrName]] <- redNucMatrix
	rm(redundantNuc)
}
for(i in 1:length(redNucList)){
	if(i == 1)	redNucMatrix <- redNucList[[i]]
	if(i != 1)	redNucMatrix <- rbind(redNucMatrix, redNucList[[i]])
}
redNucList_mm10_over500_CentralNuc <- redNucList
redNucMatrix_mm10_over500_CentralNuc <- redNucMatrix
setwd(HOMEDIR)
setwd("20211210_mouse_GRCm38_p6")
setwd("RData")
save(redNucList_mm10_over500_CentralNuc, file = "redNucList_mm10_over500_CentralNuc.RData")
save(redNucMatrix_mm10_over500_CentralNuc, file = "redNucMatrix_mm10_over500_CentralNuc.RData")


