#########################################################
## 08_prep_repH3Q85C.R
## Preparation of H3-Q85C representative nucleosomes: 20210721_nuCpos2_1.txt
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("RData")
load(file = "Chereji_wig.RData")
library(parallel)
checkMax <- function(pos, START, END){
	queryScore <- dyadsChr$score[pos]
	targetScore <- dyadsChr$score[START:END]
	if(queryScore == max(targetScore)){
		if(length(targetScore[targetScore == queryScore]) == 1) return(as.integer(1))
		if(length(targetScore[targetScore == queryScore]) > 1){
			posInStartEnd <- which((START:END) == pos)
			hitsInStartEnd <- which(targetScore == queryScore)
			if(hitsInStartEnd[1] == posInStartEnd) return(as.integer(1))
		}
		return(as.integer(0))
			
	}
	if(queryScore != max(targetScore)) return(as.integer(0))
}
countMax <- function(pos, WINDOW){
	starts <- (pos - WINDOW + 1):pos
	ends <- pos:(pos + WINDOW - 1)
	scores <- unlist(Map(f = checkMax, pos = pos, START = starts, END = ends))
	return(round(sum(scores)/length(scores), digits = 3))
}

out <- list()
WINDOW <- 107
for(i in 1:16){
	cat("i: ", i, "\n")
	dyadsChr <- wig$combined[[paste("chr", i, sep = "")]]

	colName <- paste("countMax_W", WINDOW, sep = "")
	dyadsChr[[colName]] <- as.numeric(NA)
	chrLen <- nrow(dyadsChr)
	targetPos <- WINDOW:(chrLen-WINDOW+1)
	dyadsChr[[colName]][targetPos] <- unlist(mcMap(f = countMax, 
			pos = targetPos, WINDOW = WINDOW, mc.cores = 24))
	equal1 <- nrow(subset(dyadsChr, dyadsChr[[colName]] == 1))
	cat("countMax == 1: ", equal1, "\n")
	cat("chrLen : ", chrLen, "\n")
	cat("ratio : ", equal1/chrLen, "\n")
	cat("per 147 bp : ", equal1/(chrLen/147), "\n")
	out[[paste("chr", i, sep = "")]] <- dyadsChr
	rm(dyadsChr)
} # i-loop
dyad2.Chereji <- out
rm(out)
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
dir.create("dyad_calling")
setwd("dyad_calling")
dir.create("RData_dyad")
dir.create("wiggle_dyad")
setwd("RData_dyad")
save(dyad2.Chereji, file = "dyad2.Chereji.RData")

## Generate wig files
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
wd <- getwd()
outDir1 <- "dyad_calling"
outDir2 <- "wiggle_dyad"
setwd(wd)
setwd("dyad_calling/RData_dyad")
load(file = "dyad2.Chereji.RData")
out <- dyad2.Chereji
setwd(wd)
setwd(outDir1)
setwd(outDir2)
WINDOW <- 107
prefix <- "Chereji"
w <- 1
for(j in 1:16){
	chrName <- paste("chr", j, sep = "")
	chrName2 <- as.roman(j)
	if(j == 1){
		header <- character(2)
		header[1] <- paste("track type= wiggle_0 name= \"", 
				prefix, 
				"\" description= \"dyads\"", sep = "")
		header[2] <- paste("variableStep chrom=", chrName2, " span=1", sep = "")
		write(header, file = "header.txt")
	}
	if(j > 1){
		header <- paste("variableStep chrom=", chrName2, " span=1", sep = "")
		write(header, file = "header.txt")
	}
	colName <- paste("countMax_W", WINDOW, sep = "")
	out[[j]][[colName]][is.na(out[[j]][[colName]]) == TRUE] <- 0
	write.table(out[[j]][,c(1, w+2)], file ="wig.txt", 
		sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
	file.name <- paste(prefix, "_", chrName, "_unique_dyads_w", WINDOW, ".wig", sep = "")
	command <- paste("cat header.txt wig.txt > ", file.name, sep = "")
	system(command)
	command <- "rm wig.txt header.txt"
	system(command)
}
filenames <- character(length = 16)
for(j in 1:16){
	chrName <- paste("chr", j, sep = "")
	filenames[j] <- paste(prefix, "_", chrName, "_unique_dyads_w", WINDOW, ".wig", sep = "")
}
filenames <- paste(filenames, collapse = " ", sep = "")
command <- paste("cat ", filenames, " > temp.txt", sep = "") 
system(command)
for(j in 1:16){
	chrName <- paste("chr", j, sep = "")
	command <- paste("rm ", prefix, "_", chrName, "_unique_dyads_w", WINDOW, ".wig", sep = "")
	system(command)
}
file.name <- paste(prefix, "_unique_dyads_w", WINDOW, ".wig", sep = "")
command <- paste("mv temp.txt ", file.name, sep = "") 
system(command)
file.name2 <- paste(prefix, "_unique_dyads_w", WINDOW, ".wig.tdf", sep = "")
command <- paste("igvtools toTDF ", file.name, " ", file.name2, " ",
			"/Users/hrk_kato/igv/genomes/sacCer3.genome", sep = "")
system(command)

## Generate dataframe
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
wd <- getwd()
setwd(wd)
setwd("dyad_calling/RData_dyad")
load(file = "dyad2.Chereji.RData")
out <- dyad2.Chereji
out2 <- list()
for(j in 1:16){
	chrName <- paste("chr", j, sep = "")
	temp <- out[[chrName]][out[[chrName]]$countMax_W107 == 1, ]
	out2[[chrName]] <- temp[is.na(temp$countMax_W107) == FALSE, ]
}
out3 <- data.frame()
for(j in 1:16){
	chrName <- paste("chr", j, sep = "")
	temp <- out2[[chrName]]
	temp$chr <- chrName
	temp <- temp[,c(ncol(temp), 1:(ncol(temp)-1))]
	out3 <- rbind(out3, temp)
}
out3$countMax_W107 <- NULL
Chereji_w107 <- out3
setwd(wd)
setwd("dyad_calling/RData_dyad")
save(Chereji_w107, file = "Chereji_w107.RData")


