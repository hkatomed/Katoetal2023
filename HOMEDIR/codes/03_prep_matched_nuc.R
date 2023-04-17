#########################################################
## 03_prep_matched_nuc.R
## Preparation of the matched H3-Q85C and H4-S47C datasets: 20210706_nuCpos2_1.txt
setwd(HOMEDIR)
setwd("Fuse_RData")
load(file = "nature11142_s2_sacCer3.RData")
nature11142_s2.sacCer3$term200 <- FALSE
for(i in 1:nrow(nature11142_s2.sacCer3)){
	chrName1 <- nature11142_s2.sacCer3$chr[i]
	pos <- nature11142_s2.sacCer3$pos[i]
	chrLen <- length(sc.genome.2011[[chrName1]])
	if(pos <= 200)  nature11142_s2.sacCer3$term200[i] <- TRUE
	if(pos >= (chrLen-199))  nature11142_s2.sacCer3$term200[i] <- TRUE
}
Chereji.sacCer3 <- nature11142_s2.sacCer3
Chereji.sacCer3$NCP_SN <- NULL
for(i in 1:nrow(nature11142_s2.sacCer3)){
	chrName1 <- nature11142_s2.sacCer3$chr[i]
	pos <- nature11142_s2.sacCer3$pos[i]
	posFrom <- pos - 50
	posTo <- pos + 50
	highestPos <- posFrom + which.max(wig[["combined"]][[chrName1]]$score[posFrom:posTo]) - 1
	Chereji.sacCer3$pos[i] <- highestPos
	Chereji.sacCer3$NCP_score[i] <- wig[["combined"]][[chrName1]]$score[highestPos]
}
nature11142_s2_67352 <- nature11142_s2.sacCer3[Chereji.sacCer3$NCP_score > 0 & 
				nature11142_s2.sacCer3$term200 == FALSE,]
Chereji_67352 <- Chereji.sacCer3[Chereji.sacCer3$NCP_score > 0 & 
				nature11142_s2.sacCer3$term200 == FALSE,]
nature11142_s2_67352$term200 <- NULL
Chereji_67352$term200 <- NULL
setwd(HOMEDIR)
setwd("20210706_parameters")
dir.create("RData")
setwd("RData")
save(nature11142_s2_67352, file = "nature11142_s2_67352.RData")
save(Chereji_67352, file = "Chereji_67352.RData")


