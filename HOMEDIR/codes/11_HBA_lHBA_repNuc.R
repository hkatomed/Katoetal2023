#########################################################
## 11_HBA_lHBA_repNuc.R
## Calculate HBA and local HBA of representative nucleosomes: 20210721_nuCpos2_4.txt
library("Biostrings")
setwd(HOMEDIR)
setwd("20210706_parameters")
setwd("RData")
load(file = "sc_genome_2011.RData")
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
load(file = "Chereji_w107.RData")
load(file = "Brogaard_Uniq.RData")
Brogaard_Uniq$NCP_SN <- NULL
colnames(Brogaard_Uniq)[3] <- "score"
library(Biostrings)
library(nuCpos)
library(parallel)
getSeq <- function(chrName, pos){
	left <- pos -73
	right <- pos +73
	testSeq <- sc.genome.2011[[chrName]][left:right]
	return(as.character(testSeq))
}
Chereji_w107$seq <- unlist(mcMap(f = getSeq, chrName = Chereji_w107$chr, 
					pos = Chereji_w107$pos, mc.cores = 24), 
					use.names = FALSE)
Brogaard_Uniq$seq <- unlist(mcMap(f = getSeq, chrName = Brogaard_Uniq$chr, 
					pos = Brogaard_Uniq$pos, mc.cores = 24), 
					use.names = FALSE)
Chereji_w107$HBA_sc <- unlist(mcMap(f = HBA, inseq = Chereji_w107$seq,
					species = "sc", silent = TRUE, mc.cores = 24), 
					use.names = FALSE)
Brogaard_Uniq$HBA_sc <- unlist(mcMap(f = HBA, inseq = Brogaard_Uniq$seq,
					species = "sc", silent = TRUE, mc.cores = 24), 
					use.names = FALSE)
Chereji_w107$lA_sc <- as.numeric(NA)
Chereji_w107$lB_sc <- as.numeric(NA)
Chereji_w107$lC_sc <- as.numeric(NA)
Chereji_w107$lD_sc <- as.numeric(NA)
Chereji_w107$lE_sc <- as.numeric(NA)
Chereji_w107$lF_sc <- as.numeric(NA)
Chereji_w107$lG_sc <- as.numeric(NA)
Chereji_w107$lH_sc <- as.numeric(NA)
Chereji_w107$lI_sc <- as.numeric(NA)
Chereji_w107$lJ_sc <- as.numeric(NA)
Chereji_w107$lK_sc <- as.numeric(NA)
Chereji_w107$lL_sc <- as.numeric(NA)
Chereji_w107$lM_sc <- as.numeric(NA)
Chereji_w107[, 6:18] <- matrix(unlist(
			mcMap(f = localHBA, inseq = Chereji_w107$seq,
				species = "sc", silent = TRUE, mc.cores = 24), 
				use.names = FALSE), byrow = TRUE, ncol = 13)
Brogaard_Uniq$lA_sc <- as.numeric(NA)
Brogaard_Uniq$lB_sc <- as.numeric(NA)
Brogaard_Uniq$lC_sc <- as.numeric(NA)
Brogaard_Uniq$lD_sc <- as.numeric(NA)
Brogaard_Uniq$lE_sc <- as.numeric(NA)
Brogaard_Uniq$lF_sc <- as.numeric(NA)
Brogaard_Uniq$lG_sc <- as.numeric(NA)
Brogaard_Uniq$lH_sc <- as.numeric(NA)
Brogaard_Uniq$lI_sc <- as.numeric(NA)
Brogaard_Uniq$lJ_sc <- as.numeric(NA)
Brogaard_Uniq$lK_sc <- as.numeric(NA)
Brogaard_Uniq$lL_sc <- as.numeric(NA)
Brogaard_Uniq$lM_sc <- as.numeric(NA)
Brogaard_Uniq[, 6:18] <- matrix(unlist(
			mcMap(f = localHBA, inseq = Brogaard_Uniq$seq,
				species = "sc", silent = TRUE, mc.cores = 24), 
				use.names = FALSE), byrow = TRUE, ncol = 13)
getNucName <- function(chr, pos){
	return(paste(chr, "_", pos, sep = ""))
}
Chereji_w107$name <- unlist(mcMap(f = getNucName, chr = Chereji_w107$chr, 
				pos = Chereji_w107$pos, mc.cores = 24))
Brogaard_Uniq$name <- unlist(mcMap(f = getNucName, chr = Brogaard_Uniq$chr, 
				pos = Brogaard_Uniq$pos, mc.cores = 24))
prefixs <- c("Chereji_w107", "Brogaard_Uniq")
for(p in 1:length(prefixs)){
	prefix <- prefixs[p]
	cat("# ", prefix, "\n")
	temp <- get(prefix)
	cat("#  Input:  ", nrow(temp), "nucleosomes\n")
	REMOVE1 <- round(nrow(temp)/10)
	REMOVE2 <- round(nrow(temp)/100)
	remove_low <- temp[order(temp$score), ][1:REMOVE1,]
	remove_high <- temp[order(temp$score), ][(nrow(temp)-REMOVE2+1):nrow(temp),]
	temp <- temp[order(temp$score), ][(REMOVE1+1):(nrow(temp)-REMOVE2),]
	PERGROUP <- round(nrow(temp)/9)
	cat("#  original:    ", nrow(Chereji_w107), "\n")
	cat("#  temp:        ", nrow(temp), "\n")
	cat("#  remove_high: ", nrow(remove_high), "\n")
	cat("#  remove_low:  ", nrow(remove_low), "\n")
	cat("#  PERGROUP:    ", PERGROUP, "\n")
	temp <- temp[order(temp$score, decreasing = TRUE), ]
	InObjNames <- prefix
	i <- 1
	OutObjNames <- character(length = 0)
	for(g in 1:9){
		start <- (g-1)*PERGROUP+1
		end <- g*PERGROUP
		if(g == 9) end <- nrow(temp)
		grouped <- temp[start:end, ]
		OutObjName <- paste(InObjNames[i], "_sg", g, sep = "")
		assign(OutObjName, grouped)
		OutObjNames <- c(OutObjNames, OutObjName)
		cat(paste("# sg", g, ": ", nrow(grouped), "\n", sep = ""))
	}
	temp <- get(prefix)
	temp$dyadScoreGroup <- "NA"
	for(i in 1:nrow(remove_high)){
		temp$dyadScoreGroup[temp$name == remove_high$name[i]] <- "rmHigh"
	}
	for(i in 1:nrow(remove_low)){
		temp$dyadScoreGroup[temp$name == remove_low$name[i]] <- "rmLow"
	}
	for(g in 1:9){
		groupName <- paste("sg", g, sep = "")
		cat(g, ",", sep = "")
		cat(nrow(get(OutObjNames[g])), "\n")
		for(i in 1:nrow(get(OutObjNames[g]))){
			temp$dyadScoreGroup[temp$name == get(OutObjNames[g])$name[i]] <- groupName
		}
	}
	temp$dyadScoreGroup <- factor(temp$dyadScoreGroup)
	assign(prefix, temp)
}
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
save(Chereji_w107, file = "Chereji_w107_add2.RData")
save(Brogaard_Uniq, file = "Brogaard_Uniq_add2.RData")


