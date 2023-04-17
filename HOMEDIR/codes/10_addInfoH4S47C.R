#########################################################
## 10_addInfoH4S47C.R
## Assign transcription status to the H4-S47C dataset: 20210721_nuCpos2_2.txt
setwd(HOMEDIR)
setwd("Fuse_RData")
load(file = "nature11142_s2_sacCer3.RData")
nature11142_s2.sacCer3$strand <- NULL
library("Biostrings")
setwd(HOMEDIR)
setwd("20210706_parameters")
setwd("RData")
load(file = "sc_genome_2011.RData")
chr.lengths <- width(sc.genome.2011)
chr.names <- names(sc.genome.2011)
out <- nature11142_s2.sacCer3[0,]
term <- nature11142_s2.sacCer3[0,]
for(c in 1:16){
	chr.name <- chr.names[c]
	chr.length <- chr.lengths[c]
	temp <- subset(nature11142_s2.sacCer3, chr == chr.name)
	temp1 <- subset(temp, pos > 300 & pos < (chr.length - 299))
	temp2 <- subset(temp, pos <= 300 | pos >= (chr.length - 299))
	out <- rbind(out, temp1)
	term <- rbind(term, temp2)
}
Brogaard_Uniq <- out
Brogaard_Uniq_term <- term
nrow(Brogaard_Uniq)
nrow(Brogaard_Uniq_term)
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
save(Brogaard_Uniq, file = "Brogaard_Uniq.RData")
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
load(file = "Brogaard_Uniq.RData")
setwd(HOMEDIR)
setwd("transcriptome")
setwd("RData")
load(file = "R64_2_1_GFF.RData")
load(file = "R64_2_1_sc_gff_gene_g1to13.RData")
sc.gff.gene.mRNA <- subset(R64_2_1_gene, chr != "chrmt")
Brogaard_Uniq$mRNA <- FALSE
Brogaard_Uniq$gene_id <- ""
Brogaard_Uniq$RPKB <- as.numeric(NA)
Brogaard_Uniq$strand <- ""
for(i in 1:nrow(sc.gff.gene.mRNA)){
	LEFT <- sc.gff.gene.mRNA$left[i] - 73
	RIGHT <- sc.gff.gene.mRNA$right[i] + 73
	CHR <- sc.gff.gene.mRNA$chr[i]
	RPKB <- sc.gff.gene.mRNA$RPKB[i]
	GENE_ID <- sc.gff.gene.mRNA$gene_id[i]
	STRAND <- sc.gff.gene.mRNA$strand[i]
	INDEX <- Brogaard_Uniq$chr == CHR & 
		Brogaard_Uniq$pos >= LEFT & 
		Brogaard_Uniq$pos <= RIGHT
	if(nrow(Brogaard_Uniq[INDEX,]) > 0){
		Brogaard_Uniq[INDEX,]$mRNA <- TRUE
		if(is.na(RPKB) == FALSE){	
			rpkb <- Brogaard_Uniq[INDEX,]$RPKB
			for(j in 1:length(rpkb)){
				if(is.na(rpkb[j]) == TRUE) rpkb[j] <- RPKB
				if(rpkb[j] < RPKB) rpkb[j] <- RPKB
			}
			Brogaard_Uniq[INDEX,]$RPKB <- RPKB
		}				
		gene_id <- Brogaard_Uniq[INDEX,]$gene_id
		for(j in 1:length(gene_id)){
			if(gene_id[j] == "") gene_id[j] <- GENE_ID
			if(gene_id[j] != "" & gene_id[j] != GENE_ID) gene_id[j] <- 
				paste(gene_id[j], ":", GENE_ID, sep = "", collapse = "")
		}
		Brogaard_Uniq[INDEX,]$gene_id <- gene_id
		strand <- Brogaard_Uniq[INDEX,]$strand
		for(j in 1:length(strand)){
			if(strand[j] == "") strand[j] <- STRAND
			if(strand[j] != "" & strand[j] != STRAND) strand[j] <- 
				paste(strand[j], ":", STRAND, sep = "", collapse = "")
		}
		Brogaard_Uniq[INDEX,]$strand <- strand
	}
	if(i%%100 == 0)	cat(i, ", ", sep = "")
}
Brogaard_Uniq$strand[nchar(Brogaard_Uniq$strand) > 1] <- "+/-"
Brogaard_Uniq$gRPKB <- ""
for(i in 13:1){
	group <- paste("g", i, sep = "")
	cat(group, " j:", sep = "")
	obj.name <- paste("sc.gff.gene.", group, sep = "")
	for(j in 1:nrow(get(obj.name))){
		LEFT <- get(obj.name)$left[j] - 73
		RIGHT <- get(obj.name)$right[j] + 73
		CHR <-  get(obj.name)$chr[j]
		if(nrow(Brogaard_Uniq[Brogaard_Uniq$chr == CHR & 
			Brogaard_Uniq$pos >= LEFT & 
			Brogaard_Uniq$pos <= RIGHT,]) > 0){
		Brogaard_Uniq[Brogaard_Uniq$chr == CHR & 
			Brogaard_Uniq$pos >= LEFT & 
			Brogaard_Uniq$pos <= RIGHT,]$gRPKB <- group
		}
		if(j%%100 == 0)	cat(j, ", ", sep = "")
	}
	cat("\n")
}
sc.gff.gene.mRNA$P1 <- as.integer(NA)
sc.gff.gene.mRNA$P2 <- as.integer(NA)
sc.gff.gene.mRNA$P3 <- as.integer(NA)
sc.gff.gene.mRNA$P4 <- as.integer(NA)
sc.gff.gene.mRNA$P5 <- as.integer(NA)
sc.gff.gene.mRNA$P6 <- as.integer(NA)
sc.gff.gene.mRNA$P7 <- as.integer(NA)
sc.gff.gene.mRNA$P8 <- as.integer(NA)
sc.gff.gene.mRNA$P9 <- as.integer(NA)
noHit <- 0
for(i in 1:nrow(sc.gff.gene.mRNA)){
	if(i%%100 == 0) cat(i, "/", sep = "")
	LEFT <- sc.gff.gene.mRNA$left[i] - 73
	RIGHT <- sc.gff.gene.mRNA$right[i] + 73
	CHR <- sc.gff.gene.mRNA$chr[i]
	Strand <- sc.gff.gene.mRNA$strand[i]
	if(nrow(Brogaard_Uniq[Brogaard_Uniq$chr == CHR & 
		Brogaard_Uniq$pos >= LEFT & 
		Brogaard_Uniq$pos <= RIGHT,]) > 0){
		POS <- Brogaard_Uniq[Brogaard_Uniq$chr == CHR & 
			Brogaard_Uniq$pos >= LEFT & 
			Brogaard_Uniq$pos <= RIGHT,]$pos
		if(Strand == "-") POS <- POS[length(POS):1]
		if(length(POS) >= 1) sc.gff.gene.mRNA$P1[i] <- POS[1]
		if(length(POS) >= 2) sc.gff.gene.mRNA$P2[i] <- POS[2]
		if(length(POS) >= 3) sc.gff.gene.mRNA$P3[i] <- POS[3]
		if(length(POS) >= 4) sc.gff.gene.mRNA$P4[i] <- POS[4]
		if(length(POS) >= 5) sc.gff.gene.mRNA$P5[i] <- POS[5]
		if(length(POS) >= 6) sc.gff.gene.mRNA$P6[i] <- POS[6]
		if(length(POS) >= 7) sc.gff.gene.mRNA$P7[i] <- POS[7]
		if(length(POS) >= 8) sc.gff.gene.mRNA$P8[i] <- POS[8]
		if(length(POS) >= 9) sc.gff.gene.mRNA$P9[i] <- POS[9]
	}else{
		cat("No hit: i=", i, "\n", sep = "")
		noHit <- noHit + 1
	}
}
Brogaard_Uniq$P1 <- FALSE
Brogaard_Uniq$P2 <- FALSE
Brogaard_Uniq$P3 <- FALSE
Brogaard_Uniq$P4 <- FALSE
Brogaard_Uniq$P5 <- FALSE
Brogaard_Uniq$P6 <- FALSE
Brogaard_Uniq$P7 <- FALSE
Brogaard_Uniq$P8 <- FALSE
Brogaard_Uniq$P9 <- FALSE
for(i in 1:nrow(sc.gff.gene.mRNA)){
	CHR <- sc.gff.gene.mRNA$chr[i]
	P1 <- sc.gff.gene.mRNA$P1[i]
	if(is.na(P1) == FALSE)	Brogaard_Uniq[Brogaard_Uniq$chr == CHR & 
			Brogaard_Uniq$pos == P1,]$P1 <- TRUE
	P2 <- sc.gff.gene.mRNA$P2[i]
	if(is.na(P2) == FALSE)	Brogaard_Uniq[Brogaard_Uniq$chr == CHR & 
			Brogaard_Uniq$pos == P2,]$P2 <- TRUE
	P3 <- sc.gff.gene.mRNA$P3[i]
	if(is.na(P3) == FALSE)	Brogaard_Uniq[Brogaard_Uniq$chr == CHR & 
			Brogaard_Uniq$pos == P3,]$P3 <- TRUE
	P4 <- sc.gff.gene.mRNA$P4[i]
	if(is.na(P4) == FALSE)	Brogaard_Uniq[Brogaard_Uniq$chr == CHR & 
			Brogaard_Uniq$pos == P4,]$P4 <- TRUE
	P5 <- sc.gff.gene.mRNA$P5[i]
	if(is.na(P5) == FALSE)	Brogaard_Uniq[Brogaard_Uniq$chr == CHR & 
			Brogaard_Uniq$pos == P5,]$P5 <- TRUE
	P6 <- sc.gff.gene.mRNA$P6[i]
	if(is.na(P6) == FALSE)	Brogaard_Uniq[Brogaard_Uniq$chr == CHR & 
			Brogaard_Uniq$pos == P6,]$P6 <- TRUE

	P7 <- sc.gff.gene.mRNA$P7[i]
	if(is.na(P7) == FALSE)	Brogaard_Uniq[Brogaard_Uniq$chr == CHR & 
			Brogaard_Uniq$pos == P7,]$P7 <- TRUE
	P8 <- sc.gff.gene.mRNA$P8[i]
	if(is.na(P8) == FALSE)	Brogaard_Uniq[Brogaard_Uniq$chr == CHR & 
			Brogaard_Uniq$pos == P8,]$P8 <- TRUE
	P9 <- sc.gff.gene.mRNA$P9[i]
	if(is.na(P9) == FALSE)	Brogaard_Uniq[Brogaard_Uniq$chr == CHR & 
			Brogaard_Uniq$pos == P9,]$P9 <- TRUE
}
Brogaard_Uniq$dist <- as.integer(NA)
for(i in 2:nrow(Brogaard_Uniq)){
	if(Brogaard_Uniq$chr[i-1] == Brogaard_Uniq$chr[i]){
		Brogaard_Uniq$dist[i] <- 
			Brogaard_Uniq$pos[i] - Brogaard_Uniq$pos[i-1]
	}
}
sc.gff.gene.mRNA$gRPKB <- ""
for(i in 1:13){
	group <- paste("g", i, sep = "")
	obj.name <- paste("sc.gff.gene.", group, sep = "")
	for(j in 1:nrow(get(obj.name))){
		gene_id <- get(obj.name)$gene_id[j]
		sc.gff.gene.mRNA[sc.gff.gene.mRNA$gene_id == gene_id,]$gRPKB <- group
	}
}
gRPKB <- sc.gff.gene.mRNA$gRPKB
numeric.gRPKB <- numeric(length(gRPKB))
for(i in 1:length(gRPKB)){
	if(gRPKB[i] == "")	numeric.gRPKB[i] <- 0
	if(gRPKB[i] != "")	numeric.gRPKB[i] <- as.numeric(strsplit(gRPKB[i], split = "g")[[1]][2])
}
sc.gff.gene.mRNA <- sc.gff.gene.mRNA[order(numeric.gRPKB, decreasing = TRUE),]
Brogaard_Uniq$mRNA <- FALSE
Brogaard_Uniq$gene_id <- ""
Brogaard_Uniq$RPKB <- as.numeric(NA)
Brogaard_Uniq$gRPKB <- ""
Brogaard_Uniq$strand <- ""
for(i in 1:nrow(sc.gff.gene.mRNA)){
	Strand <- sc.gff.gene.mRNA$strand[i]
	LEFT <- sc.gff.gene.mRNA$left[i] - 73
	RIGHT <- sc.gff.gene.mRNA$right[i] + 73
	P1 <- sc.gff.gene.mRNA$P1[i]
	if(Strand == "+" & is.na(P1) == FALSE)	LEFT <- P1
	if(Strand == "-" & is.na(P1) == FALSE)	RIGHT <- P1
	CHR <- sc.gff.gene.mRNA$chr[i]
	RPKB <- sc.gff.gene.mRNA$RPKB[i]
	gRPKB <- sc.gff.gene.mRNA$gRPKB[i]
	GENE_ID <- sc.gff.gene.mRNA$gene_id[i]
	INDEX <- Brogaard_Uniq$chr == CHR & 
		Brogaard_Uniq$pos >= LEFT & 
		Brogaard_Uniq$pos <= RIGHT
	if(nrow(Brogaard_Uniq[INDEX,]) > 0){
		Brogaard_Uniq[INDEX,]$mRNA <- TRUE
		rpkb <- Brogaard_Uniq[INDEX,]$RPKB
		for(j in 1:length(rpkb)){
			if(is.na(rpkb[j]) == TRUE) rpkb[j] <- RPKB
			if(is.na(RPKB) == FALSE){
				if(rpkb[j] < RPKB) rpkb[j] <- RPKB
			}
		}
		Brogaard_Uniq[INDEX,]$RPKB <- rpkb
		if(gRPKB != "")	Brogaard_Uniq[INDEX,]$gRPKB <- gRPKB
		gene_id <- Brogaard_Uniq[INDEX,]$gene_id
		for(j in 1:length(gene_id)){
			if(gene_id[j] == "") gene_id[j] <- GENE_ID
			if(gene_id[j] != "" & gene_id[j] != GENE_ID) gene_id[j] <- 
				paste(gene_id[j], ":", GENE_ID, sep = "", collapse = "")
		}
		Brogaard_Uniq[INDEX,]$gene_id <- gene_id
		strand <- Brogaard_Uniq[INDEX,]$strand
		for(j in 1:length(strand)){
			if(strand[j] == "") strand[j] <- Strand
			if(strand[j] != "" & strand[j] != Strand) strand[j] <- 
				paste(strand[j], ":", Strand, sep = "", collapse = "")
		}
		Brogaard_Uniq[INDEX,]$strand <- strand
	}
	if(i%%100 == 0)	cat(i, ", ", sep = "")
}
Brogaard_Uniq$strand[nchar(Brogaard_Uniq$strand) > 1] <- "+/-"
Brogaard_Uniq$others <- FALSE
for(i in 1:nrow(R64_2_1_others)){
	LEFT <- R64_2_1_others$left[i] - 73
	RIGHT <- R64_2_1_others$right[i] + 73
	CHR <- R64_2_1_others$chr[i]
	INDEX <- Brogaard_Uniq$chr == CHR & 
		Brogaard_Uniq$pos >= LEFT & 
		Brogaard_Uniq$pos <= RIGHT
	if(nrow(Brogaard_Uniq[INDEX,]) > 0){
		Brogaard_Uniq[INDEX,]$others <- TRUE
	}
	if(i%%100 == 0)	cat(i, ", ", sep = "")
}
Brogaard_Uniq$noMatch <- FALSE
for(i in 1:nrow(Brogaard_Uniq)){
	if(Brogaard_Uniq$mRNA[i] == FALSE & 
		Brogaard_Uniq$others[i] == FALSE){
		Brogaard_Uniq$noMatch[i] <- TRUE
	}
}
Brogaard_Uniq$type <- ""
for(i in 1:nrow(R64_2_1)){
	LEFT <- R64_2_1$left[i] - 73
	RIGHT <- R64_2_1$right[i] + 73
	CHR <- R64_2_1$chr[i]
	TYPE <- R64_2_1$type[i]
	INDEX <- Brogaard_Uniq$chr == CHR & 
		Brogaard_Uniq$pos >= LEFT & 
		Brogaard_Uniq$pos <= RIGHT
	if(nrow(Brogaard_Uniq[INDEX,]) > 0){
		type <- Brogaard_Uniq[INDEX,]$type
		for(j in 1:length(type)){
			if(type[j] == "") type[j] <- TYPE
			if(type[j] != "" & type[j] != TYPE) type[j] <- 
				paste(type[j], ":", TYPE, sep = "", collapse = "")
		}
		Brogaard_Uniq[INDEX,]$type <- type
	}
	if(i%%100 == 0)	cat(i, ", ", sep = "")
}
geneNames <- paste(Brogaard_Uniq$gene_id, collapse = ":", sep = "")
geneNames <- strsplit(geneNames, split = ":")[[1]]
geneNames <- geneNames[geneNames != ""]
countTotalNuc <- function(GENEID = "YAL068C"){
	if(length(strsplit(GENEID, split = ":")[[1]]) >= 2) return(NA)
	temp <- length(geneNames[geneNames == GENEID])
	if(temp == 0) return(NA)
	if(temp > 0) return(temp)
}
library(parallel)
Brogaard_Uniq$totalNuc <- as.integer(NA)
Brogaard_Uniq$totalNuc <- unlist(mcMap(f = countTotalNuc, GENEID = Brogaard_Uniq$gene_id, 
				mc.cores = 24), use.names = FALSE)
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
save(Brogaard_Uniq, file = "Brogaard_Uniq_add1.RData")


