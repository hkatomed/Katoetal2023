#########################################################
## 04_const_statmodels.R
## Construction of statistical models: 20210706_nuCpos2_1.txt, 20210708_nuCpos2_1.txt(33L)
setwd(HOMEDIR)
setwd("20210706_parameters")
setwd("RData")
load(file = "nature11142_s2_67352.RData")
load(file = "Chereji_67352.RData")
load(file = "sc_genome_2011.RData")
library("parallel")
library("Biostrings")
library("TTR")
nature11142_s2_67352.147 <- add.nuc.seq(dyad.table = nature11142_s2_67352, genome = sc.genome.2011, 
				species = "sc", nuc.size = 147)
Chereji_67352.147 <- add.nuc.seq(dyad.table = Chereji_67352, genome = sc.genome.2011, 
				species = "sc", nuc.size = 147)
H4S47CSeqs <- nature11142_s2_67352.147$seq
H3Q85CSeqs <- Chereji_67352.147$seq
mergedSeqs <- character(length = 67352)
for(i in 1:67352){
	leftSeq <- substr(H4S47CSeqs[i], start = 1, stop = 57)
	centerSeq <- substr(H3Q85CSeqs[i], start = 58, stop = 90)
	rightSeq <- substr(H4S47CSeqs[i], start = 91, stop = 147)
	merged33Seqs[i] <- paste(leftSeq, centerSeq, rightSeq, collapse = "", sep = "")
}
nature11142_s2_67352.147.nDSS <- append(DNAStringSet(nature11142_s2_67352.147$seq), 
			reverseComplement(DNAStringSet(nature11142_s2_67352.147$seq)))
Chereji_67352.147.nDSS <- append(DNAStringSet(Chereji_67352.147$seq), 
			reverseComplement(DNAStringSet(Chereji_67352.147$seq)))
merged33_67352.147.nDSS <- append(DNAStringSet(merged33Seqs), 
			reverseComplement(DNAStringSet(merged33Seqs)))
nature11142_s2_67352.147.freqN <- oligonucleotideFrequency(
				subseq(nature11142_s2_67352.147.nDSS, start = 1, end = 1), 
				width = 1, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE)
nature11142_s2_67352.147.tranN1 <- oligonucleotideTransitions(
				subseq(nature11142_s2_67352.147.nDSS, start = 1, end = 2), 
				left = 1, as.prob = TRUE)
nature11142_s2_67352.147.tranN2 <- oligonucleotideTransitions(
				subseq(nature11142_s2_67352.147.nDSS, start = 1, end = 3), 
				left = 2, as.prob = TRUE)
nature11142_s2_67352.147.tranN3 <- oligonucleotideTransitions(
				subseq(nature11142_s2_67352.147.nDSS, start = 1, end = 4), 
				left = 3, as.prob = TRUE)
nature11142_s2_67352.147.freqN2 <- getFreqN2(nature11142_s2_67352.147.freqN, nature11142_s2_67352.147.tranN1)
nature11142_s2_67352.147.freqN3 <- getFreqN3(nature11142_s2_67352.147.freqN2, nature11142_s2_67352.147.tranN2)
nature11142_s2_67352.147.freqN4 <- getFreqN4(nature11142_s2_67352.147.freqN3, nature11142_s2_67352.147.tranN3)
nature11142_s2_67352.147.freqN4 <- round(nature11142_s2_67352.147.freqN4, digits = 5)
nature11142_s2_67352.147.tranN4 <- getTranN4(nature11142_s2_67352.147.nDSS, nuc.size = 147)
nature11142_s2_67352.147.tranN4 <- round(nature11142_s2_67352.147.tranN4, digits = 3)
nature11142_s2_67352.147.tranN4_SMA <- getTranN4_SMA(nature11142_s2_67352.147.nDSS, nuc.size = 147)
nature11142_s2_67352.147.tranN4_SMA <- round(nature11142_s2_67352.147.tranN4_SMA, digits = 3)
Chereji_67352.147.freqN <- oligonucleotideFrequency(
				subseq(Chereji_67352.147.nDSS, start = 1, end = 1), 
				width = 1, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE)
Chereji_67352.147.tranN1 <- oligonucleotideTransitions(
				subseq(Chereji_67352.147.nDSS, start = 1, end = 2), 
				left = 1, as.prob = TRUE)
Chereji_67352.147.tranN2 <- oligonucleotideTransitions(
				subseq(Chereji_67352.147.nDSS, start = 1, end = 3), 
				left = 2, as.prob = TRUE)
Chereji_67352.147.tranN3 <- oligonucleotideTransitions(
				subseq(Chereji_67352.147.nDSS, start = 1, end = 4), 
				left = 3, as.prob = TRUE)
Chereji_67352.147.freqN2 <- getFreqN2(Chereji_67352.147.freqN, Chereji_67352.147.tranN1)
Chereji_67352.147.freqN3 <- getFreqN3(Chereji_67352.147.freqN2, Chereji_67352.147.tranN2)
Chereji_67352.147.freqN4 <- getFreqN4(Chereji_67352.147.freqN3, Chereji_67352.147.tranN3)
Chereji_67352.147.freqN4 <- round(Chereji_67352.147.freqN4, digits = 5)
Chereji_67352.147.tranN4 <- getTranN4(Chereji_67352.147.nDSS, nuc.size = 147)
Chereji_67352.147.tranN4 <- round(Chereji_67352.147.tranN4, digits = 3)
Chereji_67352.147.tranN4_SMA <- getTranN4_SMA(Chereji_67352.147.nDSS, nuc.size = 147)
Chereji_67352.147.tranN4_SMA <- round(Chereji_67352.147.tranN4_SMA, digits = 3)
merged33_67352.147.freqN <- oligonucleotideFrequency(
				subseq(merged33_67352.147.nDSS, start = 1, end = 1), 
				width = 1, as.prob = TRUE, 
				simplify.as = "collapsed", as.array = FALSE)
merged33_67352.147.tranN1 <- oligonucleotideTransitions(
				subseq(merged33_67352.147.nDSS, start = 1, end = 2), 
				left = 1, as.prob = TRUE)
merged33_67352.147.tranN2 <- oligonucleotideTransitions(
				subseq(merged33_67352.147.nDSS, start = 1, end = 3), 
				left = 2, as.prob = TRUE)
merged33_67352.147.tranN3 <- oligonucleotideTransitions(
				subseq(merged33_67352.147.nDSS, start = 1, end = 4), 
				left = 3, as.prob = TRUE)
merged33_67352.147.freqN2 <- getFreqN2(merged33_67352.147.freqN, merged33_67352.147.tranN1)
merged33_67352.147.freqN3 <- getFreqN3(merged33_67352.147.freqN2, merged33_67352.147.tranN2)
merged33_67352.147.freqN4 <- getFreqN4(merged33_67352.147.freqN3, merged33_67352.147.tranN3)
merged33_67352.147.freqN4 <- round(merged33_67352.147.freqN4, digits = 5)
merged33_67352.147.tranN4 <- getTranN4(merged33_67352.147.nDSS, nuc.size = 147)
merged33_67352.147.tranN4 <- round(merged33_67352.147.tranN4, digits = 3)
merged33_67352.147.tranN4_SMA <- getTranN4_SMA(merged33_67352.147.nDSS, nuc.size = 147)
merged33_67352.147.tranN4_SMA <- round(merged33_67352.147.tranN4_SMA, digits = 3)
save(nature11142_s2_67352.147.nDSS, file = "nature11142_s2_67352_147_nDSS.RData")
save(nature11142_s2_67352.147, file = "nature11142_s2_67352_147.RData")
save(nature11142_s2_67352.147.freqN, file = "nature11142_s2_67352_147_freqN.RData")
save(nature11142_s2_67352.147.freqN4, file = "nature11142_s2_67352_147_freqN4.RData")
save(nature11142_s2_67352.147.tranN1, file = "nature11142_s2_67352_147_tranN1.RData")
save(nature11142_s2_67352.147.tranN2, file = "nature11142_s2_67352_147_tranN2.RData")
save(nature11142_s2_67352.147.tranN3, file = "nature11142_s2_67352_147_tranN3.RData")
save(nature11142_s2_67352.147.tranN4, file = "nature11142_s2_67352_147_tranN4.RData")
save(nature11142_s2_67352.147.tranN4_SMA, file = "nature11142_s2_67352_147_tranN4_SMA.RData")
save(Chereji_67352.147.nDSS, file = "Chereji_67352_147_nDSS.RData")
save(Chereji_67352.147, file = "Chereji_67352_147.RData")
save(Chereji_67352.147.freqN, file = "Chereji_67352_147_freqN.RData")
save(Chereji_67352.147.freqN4, file = "Chereji_67352_147_freqN4.RData")
save(Chereji_67352.147.tranN1, file = "Chereji_67352_147_tranN1.RData")
save(Chereji_67352.147.tranN2, file = "Chereji_67352_147_tranN2.RData")
save(Chereji_67352.147.tranN3, file = "Chereji_67352_147_tranN3.RData")
save(Chereji_67352.147.tranN4, file = "Chereji_67352_147_tranN4.RData")
save(Chereji_67352.147.tranN4_SMA, file = "Chereji_67352_147_tranN4_SMA.RData")
save(merged33_67352.147.nDSS, file = "merged33_67352_147_nDSS.RData")
save(merged33_67352.147.freqN, file = "merged33_67352_147_freqN.RData")
save(merged33_67352.147.freqN4, file = "merged33_67352_147_freqN4.RData")
save(merged33_67352.147.tranN1, file = "merged33_67352_147_tranN1.RData")
save(merged33_67352.147.tranN2, file = "merged33_67352_147_tranN2.RData")
save(merged33_67352.147.tranN3, file = "merged33_67352_147_tranN3.RData")
save(merged33_67352.147.tranN4, file = "merged33_67352_147_tranN4.RData")
save(merged33_67352.147.tranN4_SMA, file = "merged33_67352_147_tranN4_SMA.RData")
sc.genome.2011.freqL <- oligonucleotideFrequency(sc.genome.2011, width = 1, as.prob = TRUE, 
						simplify.as = "collapsed", as.array = FALSE)
sc.genome.2011.tranL <- oligonucleotideTransitions(sc.genome.2011, left = 1, as.prob = TRUE)
sc.genome.2011.tranL2 <- oligonucleotideTransitions(sc.genome.2011, left = 2, as.prob = TRUE)
sc.genome.2011.tranL3 <- oligonucleotideTransitions(sc.genome.2011, left = 3, as.prob = TRUE)
sc.genome.2011.tranL4 <- oligonucleotideTransitions(sc.genome.2011, left = 4, as.prob = TRUE)
attr(sc.genome.2011.freqL, "names") <- NULL
attr(sc.genome.2011.tranL, "dimnames") <- NULL
attr(sc.genome.2011.tranL2, "dimnames") <- NULL
attr(sc.genome.2011.tranL3, "dimnames") <- NULL
attr(sc.genome.2011.tranL4, "dimnames") <- NULL
sc.genome.2011.freqL <- round(sc.genome.2011.freqL, digits = 3)
sc.genome.2011.tranL <- round(sc.genome.2011.tranL, digits = 3)
sc.genome.2011.tranL2 <- round(sc.genome.2011.tranL2, digits = 3)
sc.genome.2011.tranL3 <- round(sc.genome.2011.tranL3, digits = 3)
sc.genome.2011.tranL4 <- round(sc.genome.2011.tranL4, digits = 3)
save(sc.genome.2011.freqL, file = "sc_genome_2011_freqL.RData")
save(sc.genome.2011.tranL, file = "sc_genome_2011_tranL.RData")
save(sc.genome.2011.tranL2, file = "sc_genome_2011_tranL2.RData")
save(sc.genome.2011.tranL3, file = "sc_genome_2011_tranL3.RData")
save(sc.genome.2011.tranL4, file = "sc_genome_2011_tranL4.RData")


