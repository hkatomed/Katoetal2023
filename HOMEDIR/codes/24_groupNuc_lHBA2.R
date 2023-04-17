#########################################################
## 24_groupNuc_lHBA2.R
## Grouping of nucleosomes with local HBA scores (ABCKLM): 20211206_nuCpos2_2.txt
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
load(file = "Chereji_w107_add5.RData")
load(file = "Brogaard_Uniq_add5.RData")
library(parallel)
prefixs <- c("Chereji_w107", "Brogaard_Uniq")
for(p in 1:length(prefixs)){
	prefix <- prefixs[p]
	cat("# ", prefix, "\n")
	temp <- get(prefix)
	cat("#  Input:  ", nrow(temp), "nucleosomes\n")
	temp <- temp[order(temp$lA_sc), ]
	PERGROUP <- round(nrow(temp)/9)
	cat("#  original:    ", nrow( get(prefix)), "\n")
	cat("#  temp:        ", nrow(temp), "\n")
	cat("#  PERGROUP:    ", PERGROUP, "\n")
	temp <- temp[order(temp$lA_sc, decreasing = TRUE), ]
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
	temp$lA_scGroup <- "NA"
	for(g in 1:9){
		groupName <- paste("sg", g, sep = "")
		cat(g, ",", sep = "")
		cat(nrow(get(OutObjNames[g])), "\n")
		for(i in 1:nrow(get(OutObjNames[g]))){
			temp$lA_scGroup[temp$name == get(OutObjNames[g])$name[i]] <- groupName
		}
	}
	temp$lA_scGroup <- factor(temp$lA_scGroup)
	assign(prefix, temp)
}
prefixs <- c("Chereji_w107", "Brogaard_Uniq")
for(p in 1:length(prefixs)){
	prefix <- prefixs[p]
	cat("# ", prefix, "\n")
	temp <- get(prefix)
	cat("#  Input:  ", nrow(temp), "nucleosomes\n")
	temp <- temp[order(temp$lB_sc), ]
	PERGROUP <- round(nrow(temp)/9)
	cat("#  original:    ", nrow( get(prefix)), "\n")
	cat("#  temp:        ", nrow(temp), "\n")
	cat("#  PERGROUP:    ", PERGROUP, "\n")
	temp <- temp[order(temp$lB_sc, decreasing = TRUE), ]
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
	temp$lB_scGroup <- "NA"
	for(g in 1:9){
		groupName <- paste("sg", g, sep = "")
		cat(g, ",", sep = "")
		cat(nrow(get(OutObjNames[g])), "\n")
		for(i in 1:nrow(get(OutObjNames[g]))){
			temp$lB_scGroup[temp$name == get(OutObjNames[g])$name[i]] <- groupName
		}
	}
	temp$lB_scGroup <- factor(temp$lB_scGroup)
	assign(prefix, temp)
}
prefixs <- c("Chereji_w107", "Brogaard_Uniq")
for(p in 1:length(prefixs)){
	prefix <- prefixs[p]
	cat("# ", prefix, "\n")
	temp <- get(prefix)
	cat("#  Input:  ", nrow(temp), "nucleosomes\n")
	temp <- temp[order(temp$lC_sc), ]
	PERGROUP <- round(nrow(temp)/9)
	cat("#  original:    ", nrow( get(prefix)), "\n")
	cat("#  temp:        ", nrow(temp), "\n")
	cat("#  PERGROUP:    ", PERGROUP, "\n")
	temp <- temp[order(temp$lC_sc, decreasing = TRUE), ]
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
	temp$lC_scGroup <- "NA"
	for(g in 1:9){
		groupName <- paste("sg", g, sep = "")
		cat(g, ",", sep = "")
		cat(nrow(get(OutObjNames[g])), "\n")
		for(i in 1:nrow(get(OutObjNames[g]))){
			temp$lC_scGroup[temp$name == get(OutObjNames[g])$name[i]] <- groupName
		}
	}
	temp$lC_scGroup <- factor(temp$lC_scGroup)
	assign(prefix, temp)
}
prefixs <- c("Chereji_w107", "Brogaard_Uniq")
for(p in 1:length(prefixs)){
	prefix <- prefixs[p]
	cat("# ", prefix, "\n")
	temp <- get(prefix)
	cat("#  Input:  ", nrow(temp), "nucleosomes\n")
	temp <- temp[order(temp$lK_sc), ]
	PERGROUP <- round(nrow(temp)/9)
	cat("#  original:    ", nrow( get(prefix)), "\n")
	cat("#  temp:        ", nrow(temp), "\n")
	cat("#  PERGROUP:    ", PERGROUP, "\n")
	temp <- temp[order(temp$lK_sc, decreasing = TRUE), ]
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
	temp$lK_scGroup <- "NA"
	for(g in 1:9){
		groupName <- paste("sg", g, sep = "")
		cat(g, ",", sep = "")
		cat(nrow(get(OutObjNames[g])), "\n")
		for(i in 1:nrow(get(OutObjNames[g]))){
			temp$lK_scGroup[temp$name == get(OutObjNames[g])$name[i]] <- groupName
		}
	}
	temp$lK_scGroup <- factor(temp$lK_scGroup)
	assign(prefix, temp)
}
prefixs <- c("Chereji_w107", "Brogaard_Uniq")
for(p in 1:length(prefixs)){
	prefix <- prefixs[p]
	cat("# ", prefix, "\n")
	temp <- get(prefix)
	cat("#  Input:  ", nrow(temp), "nucleosomes\n")
	temp <- temp[order(temp$lL_sc), ]
	PERGROUP <- round(nrow(temp)/9)
	cat("#  original:    ", nrow( get(prefix)), "\n")
	cat("#  temp:        ", nrow(temp), "\n")
	cat("#  PERGROUP:    ", PERGROUP, "\n")
	temp <- temp[order(temp$lL_sc, decreasing = TRUE), ]
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
	temp$lL_scGroup <- "NA"
	for(g in 1:9){
		groupName <- paste("sg", g, sep = "")
		cat(g, ",", sep = "")
		cat(nrow(get(OutObjNames[g])), "\n")
		for(i in 1:nrow(get(OutObjNames[g]))){
			temp$lL_scGroup[temp$name == get(OutObjNames[g])$name[i]] <- groupName
		}
	}
	temp$lL_scGroup <- factor(temp$lL_scGroup)
	assign(prefix, temp)
}
prefixs <- c("Chereji_w107", "Brogaard_Uniq")
for(p in 1:length(prefixs)){
	prefix <- prefixs[p]
	cat("# ", prefix, "\n")
	temp <- get(prefix)
	cat("#  Input:  ", nrow(temp), "nucleosomes\n")
	temp <- temp[order(temp$lM_sc), ]
	PERGROUP <- round(nrow(temp)/9)
	cat("#  original:    ", nrow( get(prefix)), "\n")
	cat("#  temp:        ", nrow(temp), "\n")
	cat("#  PERGROUP:    ", PERGROUP, "\n")
	temp <- temp[order(temp$lM_sc, decreasing = TRUE), ]
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
	temp$lM_scGroup <- "NA"
	for(g in 1:9){
		groupName <- paste("sg", g, sep = "")
		cat(g, ",", sep = "")
		cat(nrow(get(OutObjNames[g])), "\n")
		for(i in 1:nrow(get(OutObjNames[g]))){
			temp$lM_scGroup[temp$name == get(OutObjNames[g])$name[i]] <- groupName
		}
	}
	temp$lM_scGroup <- factor(temp$lM_scGroup)
	assign(prefix, temp)
}
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
save(Chereji_w107, file = "Chereji_w107_add7.RData")
save(Brogaard_Uniq, file = "Brogaard_Uniq_add7.RData")


