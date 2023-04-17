#########################################################
## 18_groupNuc_lHBA1.R
## Grouping of nucleosomes with local HBA (FGHDEIJ): 20210721_nuCpos2_9.txt
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
load(file = "Chereji_w107_add4.RData")
load(file = "Brogaard_Uniq_add4.RData")
library(parallel)
prefixs <- c("Chereji_w107", "Brogaard_Uniq")
for(p in 1:length(prefixs)){
	prefix <- prefixs[p]
	cat("# ", prefix, "\n")
	temp <- get(prefix)
	cat("#  Input:  ", nrow(temp), "nucleosomes\n")
	temp <- temp[order(temp$lF_sc), ]
	PERGROUP <- round(nrow(temp)/9)
	cat("#  original:    ", nrow( get(prefix)), "\n")
	cat("#  temp:        ", nrow(temp), "\n")
	cat("#  PERGROUP:    ", PERGROUP, "\n")
	temp <- temp[order(temp$lF_sc, decreasing = TRUE), ]
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
	temp$lF_scGroup <- "NA"
	for(g in 1:9){
		groupName <- paste("sg", g, sep = "")
		cat(g, ",", sep = "")
		cat(nrow(get(OutObjNames[g])), "\n")
		for(i in 1:nrow(get(OutObjNames[g]))){
			temp$lF_scGroup[temp$name == get(OutObjNames[g])$name[i]] <- groupName
		}
	}
	temp$lF_scGroup <- factor(temp$lF_scGroup)
	assign(prefix, temp)
}
prefixs <- c("Chereji_w107", "Brogaard_Uniq")
for(p in 1:length(prefixs)){
	prefix <- prefixs[p]
	cat("# ", prefix, "\n")
	temp <- get(prefix)
	cat("#  Input:  ", nrow(temp), "nucleosomes\n")
	temp <- temp[order(temp$lG_sc), ]
	PERGROUP <- round(nrow(temp)/9)
	cat("#  original:    ", nrow( get(prefix)), "\n")
	cat("#  temp:        ", nrow(temp), "\n")
	cat("#  PERGROUP:    ", PERGROUP, "\n")
	temp <- temp[order(temp$lG_sc, decreasing = TRUE), ]
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
	temp$lG_scGroup <- "NA"
	for(g in 1:9){
		groupName <- paste("sg", g, sep = "")
		cat(g, ",", sep = "")
		cat(nrow(get(OutObjNames[g])), "\n")
		for(i in 1:nrow(get(OutObjNames[g]))){
			temp$lG_scGroup[temp$name == get(OutObjNames[g])$name[i]] <- groupName
		}
	}
	temp$lG_scGroup <- factor(temp$lG_scGroup)
	assign(prefix, temp)
}
prefixs <- c("Chereji_w107", "Brogaard_Uniq")
for(p in 1:length(prefixs)){
	prefix <- prefixs[p]
	cat("# ", prefix, "\n")
	temp <- get(prefix)
	cat("#  Input:  ", nrow(temp), "nucleosomes\n")
	temp <- temp[order(temp$lH_sc), ]
	PERGROUP <- round(nrow(temp)/9)
	cat("#  original:    ", nrow( get(prefix)), "\n")
	cat("#  temp:        ", nrow(temp), "\n")
	cat("#  PERGROUP:    ", PERGROUP, "\n")
	temp <- temp[order(temp$lH_sc, decreasing = TRUE), ]
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
	temp$lH_scGroup <- "NA"
	for(g in 1:9){
		groupName <- paste("sg", g, sep = "")
		cat(g, ",", sep = "")
		cat(nrow(get(OutObjNames[g])), "\n")
		for(i in 1:nrow(get(OutObjNames[g]))){
			temp$lH_scGroup[temp$name == get(OutObjNames[g])$name[i]] <- groupName
		}
	}
	temp$lH_scGroup <- factor(temp$lH_scGroup)
	assign(prefix, temp)
}
prefixs <- c("Chereji_w107", "Brogaard_Uniq")
for(p in 1:length(prefixs)){
	prefix <- prefixs[p]
	cat("# ", prefix, "\n")
	temp <- get(prefix)
	cat("#  Input:  ", nrow(temp), "nucleosomes\n")
	temp <- temp[order(temp$lD_sc), ]
	PERGROUP <- round(nrow(temp)/9)
	cat("#  original:    ", nrow( get(prefix)), "\n")
	cat("#  temp:        ", nrow(temp), "\n")
	cat("#  PERGROUP:    ", PERGROUP, "\n")
	temp <- temp[order(temp$lD_sc, decreasing = TRUE), ]
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
	temp$lD_scGroup <- "NA"
	for(g in 1:9){
		groupName <- paste("sg", g, sep = "")
		cat(g, ",", sep = "")
		cat(nrow(get(OutObjNames[g])), "\n")
		for(i in 1:nrow(get(OutObjNames[g]))){
			temp$lD_scGroup[temp$name == get(OutObjNames[g])$name[i]] <- groupName
		}
	}
	temp$lD_scGroup <- factor(temp$lD_scGroup)
	assign(prefix, temp)
}
prefixs <- c("Chereji_w107", "Brogaard_Uniq")
for(p in 1:length(prefixs)){
	prefix <- prefixs[p]
	cat("# ", prefix, "\n")
	temp <- get(prefix)
	cat("#  Input:  ", nrow(temp), "nucleosomes\n")
	temp <- temp[order(temp$lE_sc), ]
	PERGROUP <- round(nrow(temp)/9)
	cat("#  original:    ", nrow( get(prefix)), "\n")
	cat("#  temp:        ", nrow(temp), "\n")
	cat("#  PERGROUP:    ", PERGROUP, "\n")
	temp <- temp[order(temp$lE_sc, decreasing = TRUE), ]
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
	temp$lE_scGroup <- "NA"
	for(g in 1:9){
		groupName <- paste("sg", g, sep = "")
		cat(g, ",", sep = "")
		cat(nrow(get(OutObjNames[g])), "\n")
		for(i in 1:nrow(get(OutObjNames[g]))){
			temp$lE_scGroup[temp$name == get(OutObjNames[g])$name[i]] <- groupName
		}
	}
	temp$lE_scGroup <- factor(temp$lE_scGroup)
	assign(prefix, temp)
}
prefixs <- c("Chereji_w107", "Brogaard_Uniq")
for(p in 1:length(prefixs)){
	prefix <- prefixs[p]
	cat("# ", prefix, "\n")
	temp <- get(prefix)
	cat("#  Input:  ", nrow(temp), "nucleosomes\n")
	temp <- temp[order(temp$lI_sc), ]
	PERGROUP <- round(nrow(temp)/9)
	cat("#  original:    ", nrow( get(prefix)), "\n")
	cat("#  temp:        ", nrow(temp), "\n")
	cat("#  PERGROUP:    ", PERGROUP, "\n")
	temp <- temp[order(temp$lI_sc, decreasing = TRUE), ]
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
	temp$lI_scGroup <- "NA"
	for(g in 1:9){
		groupName <- paste("sg", g, sep = "")
		cat(g, ",", sep = "")
		cat(nrow(get(OutObjNames[g])), "\n")
		for(i in 1:nrow(get(OutObjNames[g]))){
			temp$lI_scGroup[temp$name == get(OutObjNames[g])$name[i]] <- groupName
		}
	}
	temp$lI_scGroup <- factor(temp$lI_scGroup)
	assign(prefix, temp)
}
prefixs <- c("Chereji_w107", "Brogaard_Uniq")
for(p in 1:length(prefixs)){
	prefix <- prefixs[p]
	cat("# ", prefix, "\n")
	temp <- get(prefix)
	cat("#  Input:  ", nrow(temp), "nucleosomes\n")
	temp <- temp[order(temp$lJ_sc), ]
	PERGROUP <- round(nrow(temp)/9)
	cat("#  original:    ", nrow( get(prefix)), "\n")
	cat("#  temp:        ", nrow(temp), "\n")
	cat("#  PERGROUP:    ", PERGROUP, "\n")
	temp <- temp[order(temp$lJ_sc, decreasing = TRUE), ]
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
	temp$lJ_scGroup <- "NA"
	for(g in 1:9){
		groupName <- paste("sg", g, sep = "")
		cat(g, ",", sep = "")
		cat(nrow(get(OutObjNames[g])), "\n")
		for(i in 1:nrow(get(OutObjNames[g]))){
			temp$lJ_scGroup[temp$name == get(OutObjNames[g])$name[i]] <- groupName
		}
	}
	temp$lJ_scGroup <- factor(temp$lJ_scGroup)
	assign(prefix, temp)
}
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
save(Chereji_w107, file = "Chereji_w107_add5.RData")
save(Brogaard_Uniq, file = "Brogaard_Uniq_add5.RData")


