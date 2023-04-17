#########################################################
## 36_calcHBAexonsMouse4.R
## Calculation of HBA along protein-coding genes (mouse)
## Preparation of exon: 20211210_nuCpos2_1.txt			matchedExon.RData
## Download mm10 GTF from NCBI (GRCm38_p6) and store in
## HOMEDIR/20211210_mouse_GRCm38_p6/annotation/genome_assemblies_genome_gtf, and untar it.
setwd(HOMEDIR)
setwd("20211210_mouse_GRCm38_p6/annotation/NCBI_GTF")
setwd("genome_assemblies_genome_gtf/ncbi-genomes-2021-12-09")
GTF <- read.table(file = "GCF_000001635.26_GRCm38.p6_genomic.gtf.gz", skip = 4, sep = "\t", header = FALSE)
RefSeqSelect <- GTF[grep(pattern = "RefSeq Select", GTF$V9),]
getGeneID <- function(V9){
	GeneID <- strsplit(strsplit(V9, split = "; ")[[1]][3], split = " ")[[1]][2]
	return(GeneID)
}
getTranscriptID <- function(V9){
	TranscriptID <- strsplit(strsplit(V9, split = "; ")[[1]][2], split = " ")[[1]][2]
	return(TranscriptID)
}
library(parallel)
RefSeqSelect$GeneID <- unlist(mcMap(f = getGeneID, 
				V9 = RefSeqSelect$V9, mc.cores = 24), use.names = FALSE)
RefSeqSelect$TranscriptID <- unlist(mcMap(f = getTranscriptID, 
				V9 = RefSeqSelect$V9, mc.cores = 24), use.names = FALSE)
RefSeqSelectExon <- subset(RefSeqSelect, RefSeqSelect$V3 == "exon")
RefSeqSelectExon$GeneID <- unlist(mcMap(f = getGeneID, 
				V9 = RefSeqSelectExon$V9, mc.cores = 24), use.names = FALSE)
RefSeqSelectExon$TranscriptID <- unlist(mcMap(f = getTranscriptID, 
				V9 = RefSeqSelectExon$V9, mc.cores = 24), use.names = FALSE)
RefSeqSelectCDS <- subset(RefSeqSelect, RefSeqSelect$V3 == "CDS")
RefSeqSelectCDS$GeneID <- unlist(mcMap(f = getGeneID, 
				V9 = RefSeqSelectCDS$V9, mc.cores = 24), use.names = FALSE)
RefSeqSelectCDS$TranscriptID <- unlist(mcMap(f = getTranscriptID, 
				V9 = RefSeqSelectCDS$V9, mc.cores = 24), use.names = FALSE)
exon.GeneID <- unique(RefSeqSelectExon$GeneID)
exon.GeneID <- exon.GeneID[order(exon.GeneID)]
exon.TranscriptID <- unique(RefSeqSelectExon$TranscriptID)
exon.TranscriptID <- exon.TranscriptID[order(exon.TranscriptID)]
CDS.GeneID <- unique(RefSeqSelectCDS$GeneID)
CDS.GeneID <- CDS.GeneID[order(CDS.GeneID)]
CDS.TranscriptID <- unique(RefSeqSelectCDS$TranscriptID)
CDS.TranscriptID <- CDS.TranscriptID[order(CDS.TranscriptID)]
total.TranscriptID <- unique(c(exon.TranscriptID, CDS.TranscriptID))
total.TranscriptID <- total.TranscriptID[order(total.TranscriptID)]
compareTiD <- data.frame(TranscriptID = exon.TranscriptID, CDS = FALSE)
for(n in 1:length(CDS.TranscriptID)){
	compareTiD$CDS[compareTiD$TranscriptID == CDS.TranscriptID[n]] <- TRUE
}
getGeneName <- function(V9){
	GeneName <- strsplit(strsplit(V9, split = "; ")[[1]][1], split = " ")[[1]][2]
	return(GeneName)
}
getExonNumber_temp <- function(V9){
	temp <- strsplit(V9, split = "; ")[[1]]
	ExonNumber_temp <- temp[grep(pattern = "exon_number", temp)]
	return(ExonNumber_temp)
}
RefSeqSelectExon$GeneName <- unlist(mcMap(f = getGeneName, 
				V9 = RefSeqSelectExon$V9, mc.cores = 24), use.names = FALSE)
RefSeqSelectExon$ExonNumber_temp <- unlist(mcMap(f = getExonNumber_temp, 
				V9 = RefSeqSelectExon$V9, mc.cores = 24), use.names = FALSE)
RefSeqSelectExon$ExonNumber_temp <- NULL
getExonNumber <- function(V9){
	temp <- strsplit(V9, split = "; ")[[1]]
	temp2 <- temp[grep(pattern = "exon_number", temp)]
	ExonNumber <- as.integer(strsplit(temp2, split = " ")[[1]][2])
	return(ExonNumber)
}
RefSeqSelectExon$ExonNumber <- unlist(mcMap(f = getExonNumber, 
				V9 = RefSeqSelectExon$V9, mc.cores = 24), use.names = FALSE)
RefSeqSelectCDS$GeneName <- unlist(mcMap(f = getGeneName, 
				V9 = RefSeqSelectCDS$V9, mc.cores = 24), use.names = FALSE)
RefSeqSelectCDS$ExonNumber <- unlist(mcMap(f = getExonNumber, 
				V9 = RefSeqSelectCDS$V9, mc.cores = 24), use.names = FALSE)
RefSeqSelectExon$ExonNumber2 <- RefSeqSelectExon$ExonNumber
for(i in 1:length(exon.TranscriptID)){
	if(i%%500 == 0)	cat(i, ",", sep = "")
	TiD <- exon.TranscriptID[i]
	index <- grep(pattern = TiD, RefSeqSelectExon$TranscriptID)
	strand <- unique(RefSeqSelectExon$V7[index])
	if(length(strand) != 1)	stop("Multiple strands: ", i, sep = "")
	if(strand[1] == "-"){
		originalEN <- RefSeqSelectExon$ExonNumber[index]
		firstEN <- min(originalEN)
		lastEN <- max(originalEN)
		modifiedEN <- lastEN:firstEN
		originalEN_Rvc <- paste(originalEN[length(originalEN):1], collapse = "", sep = "")
		modifiedENc <- paste(modifiedEN, collapse = "", sep = "")
		if(originalEN_Rvc != modifiedENc){
			stop("modifiedEN not equeal: ", i, sep = "")
		}
		if(originalEN_Rvc == modifiedENc){
			RefSeqSelectExon$ExonNumber2[index] <- modifiedEN
		}
	}
}	
RefSeqSelectCDS$ExonNumber2 <- RefSeqSelectCDS$ExonNumber
for(i in c(1:10570, 10572:length(CDS.TranscriptID))){
	if(i%%500 == 0)	cat(i, ",", sep = "")
	TiD <- CDS.TranscriptID[i]
	index <- grep(pattern = TiD, RefSeqSelectCDS$TranscriptID)
	strand <- unique(RefSeqSelectCDS$V7[index])
	if(length(strand) != 1)	stop("Multiple strands: ", i, sep = "")
	if(strand[1] == "-"){
		originalEN <- RefSeqSelectCDS$ExonNumber[index]
		firstEN <- min(originalEN)
		lastEN <- max(originalEN)
		modifiedEN <- lastEN:firstEN
		originalEN_Rvc <- paste(originalEN[length(originalEN):1], collapse = "", sep = "")
		modifiedENc <- paste(modifiedEN, collapse = "", sep = "")
		if(originalEN_Rvc != modifiedENc){
			stop("modifiedEN not equeal: ", i, sep = "")
		}
		if(originalEN_Rvc == modifiedENc){
			RefSeqSelectCDS$ExonNumber2[index] <- modifiedEN
		}
	}
}
RefSeqSelectCDS$ExonNumber2 <- NULL	
RefSeqSelectExon$CDSmatch <- FALSE
for(i in 1:nrow(RefSeqSelectCDS)){
	CDS.TranscriptID <- RefSeqSelectCDS$TranscriptID[i]
	CDS.ExonNumber <- RefSeqSelectCDS$ExonNumber[i]
	
	index <- RefSeqSelectExon$TranscriptID == CDS.TranscriptID &
			RefSeqSelectExon$ExonNumber == CDS.ExonNumber
	index <- which(index == TRUE)
	if(length(index) == 1){
		Exon.chr <- RefSeqSelectExon$V1[index]
		Exon.start <- RefSeqSelectExon$V4[index]
		Exon.end <- RefSeqSelectExon$V5[index]
		CDS.chr <- RefSeqSelectCDS$V1[i]
		CDS.start <- RefSeqSelectCDS$V4[i]
		CDS.end <- RefSeqSelectCDS$V5[i]
		if(Exon.chr == CDS.chr & Exon.start == CDS.start & 
				Exon.end == CDS.end){
			RefSeqSelectExon$CDSmatch[index] <- TRUE
		}
	}
	if(length(index) != 1){
		cat("i: ", i, "index: ", index, " TiD: ", CDS.TranscriptID, "\n")
	}
}
matchedExon <- subset(RefSeqSelectExon, CDSmatch == TRUE)[, c(1, 4, 5, 7, 10:15)]
colnames(matchedExon)[c(1:4)] <- c("chrName", "left", "right", "strand")
matchedExon$length <- matchedExon$right - matchedExon$left + 1
setwd(HOMEDIR)
setwd("20211210_mouse_GRCm38_p6/NCBI_FASTA")
setwd("genome_assemblies_asm_struct")
setwd("ncbi-genomes-2021-12-10")
assemblyReport <- read.table(file = "GCF_000001635.26_GRCm38.p6_assembly_report.txt", skip = 43, header = FALSE)[, c(7, 10)]
mm10names <- matchedExon$chrName
for(n in 1:length(mm10names)){
	NCBIname <- strsplit(mm10names[n], split = " ")[[1]][1]
	mm10names[n] <- assemblyReport$V10[assemblyReport$V7 == NCBIname]
}
matchedExon$chrNameUCSC <- mm10names
matchedExon_old <- matchedExon
matchedExon <- subset(matchedExon, chrNameUCSC != "na")
getUniqueName <- function(TranscriptID, ExonNumber){
	temp <- paste(TranscriptID, "_", ExonNumber, sep = "")
	return(temp)
}
matchedExon$UniqueName <- unlist(mcMap(f = getUniqueName, 
	TranscriptID = matchedExon$TranscriptID, 
	ExonNumber = matchedExon$ExonNumber2, mc.cores = 24), use.names = FALSE)
setwd(HOMEDIR)
setwd("20211210_mouse_GRCm38_p6")
dir.create("RData")
setwd("RData")
save(RefSeqSelectExon, file = "RefSeqSelectExon.RData")
save(RefSeqSelectCDS, file = "RefSeqSelectCDS.RData")
save(exon.TranscriptID, file = "exon_TranscriptID.RData")
save(CDS.TranscriptID, file = "CDS_TranscriptID.RData")
save(matchedExon_old, file = "matchedExon_old.RData")
save(matchedExon, file = "matchedExon.RData")


