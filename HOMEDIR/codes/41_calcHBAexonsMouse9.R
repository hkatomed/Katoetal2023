#########################################################
## 41_calcHBAexonsMouse9.R
## Calculation of HBA along protein-coding genes (mouse)
## Preparation of exons: 20211217_nuCpos2_2.txt		nucSeq_over500_CentralNuc.RData
library(Biostrings)
library(parallel)
setwd(HOMEDIR)
setwd("20211210_mouse_GRCm38_p6")
setwd("RData")
load(file = "mm10.RData")
setwd(HOMEDIR)
setwd("20211210_mouse_GRCm38_p6")
setwd("RData")
load(file = "matchedExon_over500_CentralNuc.RData")
getSeqCentralNuc <- function(chrName, centralNuc, strand){
	left <- centralNuc - 500 -73
	right <- centralNuc + 500 +73
	
	SEQ <- mm10[[chrName]][left:right]
	if(strand == "-")	SEQ <- reverseComplement(SEQ)
	return(as.character(SEQ))
}
nucSeq_over500_CentralNuc <- unlist(mcMap(f = getSeqCentralNuc, 
			chrName = matchedExon_over500_CentralNuc$chrNameUCSC, 
			centralNuc = matchedExon_over500_CentralNuc$centralNuc, 
			strand = matchedExon_over500_CentralNuc$strand,
			mc.cores = 24), use.names = FALSE)
names(nucSeq_over500_CentralNuc) <- matchedExon_over500_CentralNuc$UniqueName
nucSeq_over500_CentralNuc <- DNAStringSet(nucSeq_over500_CentralNuc)
check.noN.inDSS <- function(SEQ){
	Nnum <- length(matchPattern(pattern = "N", subject = SEQ))
	if(Nnum == 0) return(TRUE)
	if(Nnum != 0) return(FALSE)
}
noN.index <- unlist(mcMap(f = check.noN.inDSS, 
			as.character(nucSeq_over500_CentralNuc), 
			mc.cores = 24), use.names = FALSE)
temp <- vmatchPattern(pattern = "N", subject = nucSeq_over500_CentralNuc)
for(t in 1:length(temp)){
	if(length(temp[[t]]) != 0){
		cat("# t:", t, " ID:", names(temp)[t], "\n", sep = "")
	}
}
setwd(HOMEDIR)
setwd("20211210_mouse_GRCm38_p6")
setwd("RData")
save(nucSeq_over500_CentralNuc, file = "nucSeq_over500_CentralNuc.RData")	


