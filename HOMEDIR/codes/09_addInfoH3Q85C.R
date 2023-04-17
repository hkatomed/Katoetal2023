#########################################################
## 09_addInfoH3Q85C.R
## Assign transcription status to the H3-Q85C dataset: 20210721_nuCpos2_2.txt
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
load(file = "Chereji_w107.RData")
setwd(HOMEDIR)
setwd("transcriptome")
setwd("RData")
load(file = "R64_2_1_GFF.RData")
load(file = "R64_2_1_sc_gff_gene_g1to13.RData")
sc.gff.gene.mRNA <- subset(R64_2_1_gene, chr != "chrmt")
Chereji_w107$mRNA <- FALSE
Chereji_w107$gene_id <- ""
Chereji_w107$RPKB <- as.numeric(NA)
Chereji_w107$strand <- ""
for(i in 1:nrow(sc.gff.gene.mRNA)){
	LEFT <- sc.gff.gene.mRNA$left[i] - 73
	RIGHT <- sc.gff.gene.mRNA$right[i] + 73
	CHR <- sc.gff.gene.mRNA$chr[i]
	RPKB <- sc.gff.gene.mRNA$RPKB[i]
	GENE_ID <- sc.gff.gene.mRNA$gene_id[i]
	STRAND <- sc.gff.gene.mRNA$strand[i]
	INDEX <- Chereji_w107$chr == CHR & 
		Chereji_w107$pos >= LEFT & 
		Chereji_w107$pos <= RIGHT
	if(nrow(Chereji_w107[INDEX,]) > 0){
		Chereji_w107[INDEX,]$mRNA <- TRUE
		if(is.na(RPKB) == FALSE){	
			rpkb <- Chereji_w107[INDEX,]$RPKB
			for(j in 1:length(rpkb)){
				if(is.na(rpkb[j]) == TRUE) rpkb[j] <- RPKB
				if(rpkb[j] < RPKB) rpkb[j] <- RPKB
			}
			Chereji_w107[INDEX,]$RPKB <- RPKB
		}				
		gene_id <- Chereji_w107[INDEX,]$gene_id
		for(j in 1:length(gene_id)){
			if(gene_id[j] == "") gene_id[j] <- GENE_ID
			if(gene_id[j] != "" & gene_id[j] != GENE_ID) gene_id[j] <- 
				paste(gene_id[j], ":", GENE_ID, sep = "", collapse = "")
		}
		Chereji_w107[INDEX,]$gene_id <- gene_id
		strand <- Chereji_w107[INDEX,]$strand
		for(j in 1:length(strand)){
			if(strand[j] == "") strand[j] <- STRAND
			if(strand[j] != "" & strand[j] != STRAND) strand[j] <- 
				paste(strand[j], ":", STRAND, sep = "", collapse = "")
		}
		Chereji_w107[INDEX,]$strand <- strand
	}
	if(i%%100 == 0)	cat(i, ", ", sep = "")
}
Chereji_w107$strand[nchar(Chereji_w107$strand) > 1] <- "+/-"
Chereji_w107$gRPKB <- ""
for(i in 13:1){
	group <- paste("g", i, sep = "")
	cat(group, " j:", sep = "")
	obj.name <- paste("sc.gff.gene.", group, sep = "")
	for(j in 1:nrow(get(obj.name))){
		LEFT <- get(obj.name)$left[j] - 73
		RIGHT <- get(obj.name)$right[j] + 73
		CHR <-  get(obj.name)$chr[j]
		if(nrow(Chereji_w107[Chereji_w107$chr == CHR & 
			Chereji_w107$pos >= LEFT & 
			Chereji_w107$pos <= RIGHT,]) > 0){
		Chereji_w107[Chereji_w107$chr == CHR & 
			Chereji_w107$pos >= LEFT & 
			Chereji_w107$pos <= RIGHT,]$gRPKB <- group
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
	if(nrow(Chereji_w107[Chereji_w107$chr == CHR & 
		Chereji_w107$pos >= LEFT & 
		Chereji_w107$pos <= RIGHT,]) > 0){
		POS <- Chereji_w107[Chereji_w107$chr == CHR & 
			Chereji_w107$pos >= LEFT & 
			Chereji_w107$pos <= RIGHT,]$pos
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
Chereji_w107$P1 <- FALSE
Chereji_w107$P2 <- FALSE
Chereji_w107$P3 <- FALSE
Chereji_w107$P4 <- FALSE
Chereji_w107$P5 <- FALSE
Chereji_w107$P6 <- FALSE
Chereji_w107$P7 <- FALSE
Chereji_w107$P8 <- FALSE
Chereji_w107$P9 <- FALSE
for(i in 1:nrow(sc.gff.gene.mRNA)){
	CHR <- sc.gff.gene.mRNA$chr[i]
	P1 <- sc.gff.gene.mRNA$P1[i]
	if(is.na(P1) == FALSE)	Chereji_w107[Chereji_w107$chr == CHR & 
			Chereji_w107$pos == P1,]$P1 <- TRUE
	P2 <- sc.gff.gene.mRNA$P2[i]
	if(is.na(P2) == FALSE)	Chereji_w107[Chereji_w107$chr == CHR & 
			Chereji_w107$pos == P2,]$P2 <- TRUE
	P3 <- sc.gff.gene.mRNA$P3[i]
	if(is.na(P3) == FALSE)	Chereji_w107[Chereji_w107$chr == CHR & 
			Chereji_w107$pos == P3,]$P3 <- TRUE
	P4 <- sc.gff.gene.mRNA$P4[i]
	if(is.na(P4) == FALSE)	Chereji_w107[Chereji_w107$chr == CHR & 
			Chereji_w107$pos == P4,]$P4 <- TRUE
	P5 <- sc.gff.gene.mRNA$P5[i]
	if(is.na(P5) == FALSE)	Chereji_w107[Chereji_w107$chr == CHR & 
			Chereji_w107$pos == P5,]$P5 <- TRUE
	P6 <- sc.gff.gene.mRNA$P6[i]
	if(is.na(P6) == FALSE)	Chereji_w107[Chereji_w107$chr == CHR & 
			Chereji_w107$pos == P6,]$P6 <- TRUE

	P7 <- sc.gff.gene.mRNA$P7[i]
	if(is.na(P7) == FALSE)	Chereji_w107[Chereji_w107$chr == CHR & 
			Chereji_w107$pos == P7,]$P7 <- TRUE
	P8 <- sc.gff.gene.mRNA$P8[i]
	if(is.na(P8) == FALSE)	Chereji_w107[Chereji_w107$chr == CHR & 
			Chereji_w107$pos == P8,]$P8 <- TRUE
	P9 <- sc.gff.gene.mRNA$P9[i]
	if(is.na(P9) == FALSE)	Chereji_w107[Chereji_w107$chr == CHR & 
			Chereji_w107$pos == P9,]$P9 <- TRUE
}
Chereji_w107$dist <- as.integer(NA)
for(i in 2:nrow(Chereji_w107)){
	if(Chereji_w107$chr[i-1] == Chereji_w107$chr[i]){
		Chereji_w107$dist[i] <- 
			Chereji_w107$pos[i] - Chereji_w107$pos[i-1]
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
Chereji_w107$mRNA <- FALSE
Chereji_w107$gene_id <- ""
Chereji_w107$RPKB <- as.numeric(NA)
Chereji_w107$gRPKB <- ""
Chereji_w107$strand <- ""
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
	INDEX <- Chereji_w107$chr == CHR & 
		Chereji_w107$pos >= LEFT & 
		Chereji_w107$pos <= RIGHT
	if(nrow(Chereji_w107[INDEX,]) > 0){
		Chereji_w107[INDEX,]$mRNA <- TRUE
		rpkb <- Chereji_w107[INDEX,]$RPKB
		for(j in 1:length(rpkb)){
			if(is.na(rpkb[j]) == TRUE) rpkb[j] <- RPKB
			if(is.na(RPKB) == FALSE){
				if(rpkb[j] < RPKB) rpkb[j] <- RPKB
			}
		}
		Chereji_w107[INDEX,]$RPKB <- rpkb
		if(gRPKB != "")	Chereji_w107[INDEX,]$gRPKB <- gRPKB
		gene_id <- Chereji_w107[INDEX,]$gene_id
		for(j in 1:length(gene_id)){
			if(gene_id[j] == "") gene_id[j] <- GENE_ID
			if(gene_id[j] != "" & gene_id[j] != GENE_ID) gene_id[j] <- 
				paste(gene_id[j], ":", GENE_ID, sep = "", collapse = "")
		}
		Chereji_w107[INDEX,]$gene_id <- gene_id
		strand <- Chereji_w107[INDEX,]$strand
		for(j in 1:length(strand)){
			if(strand[j] == "") strand[j] <- Strand
			if(strand[j] != "" & strand[j] != Strand) strand[j] <- 
				paste(strand[j], ":", Strand, sep = "", collapse = "")
		}
		Chereji_w107[INDEX,]$strand <- strand
	}
	if(i%%100 == 0)	cat(i, ", ", sep = "")
}
Chereji_w107$strand[nchar(Chereji_w107$strand) > 1] <- "+/-"
Chereji_w107$others <- FALSE
for(i in 1:nrow(R64_2_1_others)){
	LEFT <- R64_2_1_others$left[i] - 73
	RIGHT <- R64_2_1_others$right[i] + 73
	CHR <- R64_2_1_others$chr[i]
	INDEX <- Chereji_w107$chr == CHR & 
		Chereji_w107$pos >= LEFT & 
		Chereji_w107$pos <= RIGHT
	if(nrow(Chereji_w107[INDEX,]) > 0){
		Chereji_w107[INDEX,]$others <- TRUE
	}
	if(i%%100 == 0)	cat(i, ", ", sep = "")
}
Chereji_w107$noMatch <- FALSE
for(i in 1:nrow(Chereji_w107)){
	if(Chereji_w107$mRNA[i] == FALSE & 
		Chereji_w107$others[i] == FALSE){
		Chereji_w107$noMatch[i] <- TRUE
	}
}
Chereji_w107$type <- ""
for(i in 1:nrow(R64_2_1)){
	LEFT <- R64_2_1$left[i] - 73
	RIGHT <- R64_2_1$right[i] + 73
	CHR <- R64_2_1$chr[i]
	TYPE <- R64_2_1$type[i]
	INDEX <- Chereji_w107$chr == CHR & 
		Chereji_w107$pos >= LEFT & 
		Chereji_w107$pos <= RIGHT
	if(nrow(Chereji_w107[INDEX,]) > 0){
		type <- Chereji_w107[INDEX,]$type
		for(j in 1:length(type)){
			if(type[j] == "") type[j] <- TYPE
			if(type[j] != "" & type[j] != TYPE) type[j] <- 
				paste(type[j], ":", TYPE, sep = "", collapse = "")
		}
		Chereji_w107[INDEX,]$type <- type
	}
	if(i%%100 == 0)	cat(i, ", ", sep = "")
}
geneNames <- paste(Chereji_w107$gene_id, collapse = ":", sep = "")
geneNames <- strsplit(geneNames, split = ":")[[1]]
geneNames <- geneNames[geneNames != ""]
countTotalNuc <- function(GENEID = "YAL068C"){
	if(length(strsplit(GENEID, split = ":")[[1]]) >= 2) return(NA)
	temp <- length(geneNames[geneNames == GENEID])
	if(temp == 0) return(NA)
	if(temp > 0) return(temp)
}
library(parallel)
Chereji_w107$totalNuc <- as.integer(NA)
Chereji_w107$totalNuc <- unlist(mcMap(f = countTotalNuc, GENEID = Chereji_w107$gene_id, 
				mc.cores = 24), use.names = FALSE)
setwd(HOMEDIR)
setwd("20210706_Chereji_GSE97290")
setwd("dyad_calling/RData_dyad")
save(Chereji_w107, file = "Chereji_w107_add1.RData")


