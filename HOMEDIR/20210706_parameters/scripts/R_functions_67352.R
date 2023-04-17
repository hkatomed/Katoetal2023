
add.nuc.seq <- function(dyad.table, genome, species = "sc", nuc.size = 147){
	require("Biostrings")
	require("parallel")
	LEN <- (nuc.size - 1)/2
	dyad.table$seq <- as.character(NA)
	if(species == "sc"){
		dyad.table.chr1 <- subset(dyad.table, chr == "chr1" & pos > LEN & pos < (230208-LEN))
		dyad.table.chr2 <- subset(dyad.table, chr == "chr2" & pos > LEN & pos < (813178-LEN))
		dyad.table.chr3 <- subset(dyad.table, chr == "chr3" & pos > LEN & pos < (316617-LEN))
		dyad.table.chr4 <- subset(dyad.table, chr == "chr4" & pos > LEN & pos < (1531919-LEN))
		dyad.table.chr5 <- subset(dyad.table, chr == "chr5" & pos > LEN & pos < (576869-LEN))
		dyad.table.chr6 <- subset(dyad.table, chr == "chr6" & pos > LEN & pos < (270148-LEN))
		dyad.table.chr7 <- subset(dyad.table, chr == "chr7" & pos > LEN & pos < (1090947-LEN))
		dyad.table.chr8 <- subset(dyad.table, chr == "chr8" & pos > LEN & pos < (562643-LEN))
		dyad.table.chr9 <- subset(dyad.table, chr == "chr9" & pos > LEN & pos < (439885-LEN))
		dyad.table.chr10 <- subset(dyad.table, chr == "chr10" & pos > LEN & pos < (745742-LEN))
		dyad.table.chr11 <- subset(dyad.table, chr == "chr11" & pos > LEN & pos < (666454-LEN))
		dyad.table.chr12 <- subset(dyad.table, chr == "chr12" & pos > LEN & pos < (1078175-LEN))
		dyad.table.chr13 <- subset(dyad.table, chr == "chr13" & pos > LEN & pos < (924429-LEN))
		dyad.table.chr14 <- subset(dyad.table, chr == "chr14" & pos > LEN & pos < (784333-LEN))
		dyad.table.chr15 <- subset(dyad.table, chr == "chr15" & pos > LEN & pos < (1091289-LEN))
		dyad.table.chr16 <- subset(dyad.table, chr == "chr16" & pos > LEN & pos < (948062-LEN))
		dyad.table <- rbind(dyad.table.chr1, dyad.table.chr2, dyad.table.chr3, 
			dyad.table.chr4, dyad.table.chr5, dyad.table.chr6, dyad.table.chr7, 
			dyad.table.chr8, dyad.table.chr9, dyad.table.chr10, dyad.table.chr11, 
			dyad.table.chr12, dyad.table.chr13, dyad.table.chr14, dyad.table.chr15, 
			dyad.table.chr16)
	}else if(species == "sp"){
		dyad.table.chr1 <- subset(dyad.table, chr == "chr1" & pos > LEN & pos < (5579133-LEN))
		dyad.table.chr2 <- subset(dyad.table, chr == "chr2" & pos > LEN & pos < (4539804-LEN))
		dyad.table.chr3 <- subset(dyad.table, chr == "chr3" & pos > LEN & pos < (2452883-LEN))
		dyad.table <- rbind(dyad.table.chr1, dyad.table.chr2, dyad.table.chr3)
	}
	
		get.seq <- function(chromosome, pos, strand){
			left <- pos - LEN
			right <- pos + LEN
			if(strand == "+")	out.seq <- 
				as.character(genome[[chromosome]][left:right])
			if(strand == "-")	out.seq <- 
				as.character(reverseComplement(genome[[chromosome]][left:right]))
			return(out.seq)
		}

	dyad.table$seq <- unlist(mcMap(get.seq, chromosome = dyad.table$chr, pos = dyad.table$pos, 
				strand = dyad.table$strand, mc.cores = 20))
	return(dyad.table)
}

getFreqN2 <- function(freqN, tranN1){
	for(i in 1:4){
		tranN1[i,] <- tranN1[i,]*freqN[i]
	}
	return(tranN1)
}

getFreqN3 <- function(freqN2, tranN2){
	for(i in 1:nrow(tranN2)){
		tranN2[i,] <- tranN2[i,]*t(freqN2)[i]
	}
	return(tranN2)
}

getFreqN4 <- function(freqN3, tranN3){
	for(i in 1:nrow(tranN3)){
		tranN3[i,] <- tranN3[i,]*t(freqN3)[i]
	}
	return(tranN3)
}

getTranN4 <- function(nDSS, nuc.size = 147){
	require(Biostrings)
	TranN4 <- matrix(numeric(0), ncol = 4)
	for(i in 1:(nuc.size-4)){
		iTranN4 <- oligonucleotideTransitions(subseq(nDSS, start = i, end = i + 4),
				 left = 4, as.prob = TRUE)
		TranN4 <- rbind(TranN4, iTranN4)
	}
	return(TranN4)
}

getTranN4_SMA <- function(nDSS, nuc.size = 147){
	require(Biostrings)
	require(TTR)
	TranN4 <- matrix(numeric(0), nrow = 0, ncol = 4)
	TranN4.table <- matrix(numeric((nuc.size-4)*256), nrow = (nuc.size-4), ncol = 1024)
	TranN4.table.SMA <- matrix(numeric(0), nrow = (nuc.size-4), ncol = 1024)
	for(i in 1:(nuc.size-4)){
		iTranN4 <- oligonucleotideTransitions(subseq(nDSS, start = i, end = i + 4),
				 left = 4, as.prob = TRUE)
		TranN4.table[i,] <- as.numeric(t(iTranN4))
	}
	for(i in 1:1024){
		TranN4.table.SMA[,i] <- SMA(TranN4.table[,i], n = 3)
	}
	
	TranN4.table.SMA[2:(nuc.size-5),] <- TranN4.table.SMA[3:(nuc.size-4),]
	TranN4.table.SMA[1,] <- TranN4.table[1,]
	TranN4.table.SMA[(nuc.size-4),] <- TranN4.table[(nuc.size-4),]

	for(i in 1:(nuc.size-4)){
		TranN4 <- rbind(TranN4, matrix(TranN4.table.SMA[i,], byrow = TRUE, ncol = 4))
	}
	
	# SMA後の補正
	for(i in 1:nrow(TranN4)){
		TranN4[i,] <- TranN4[i,]/sum(TranN4[i,])
	}

	return(TranN4)
}



