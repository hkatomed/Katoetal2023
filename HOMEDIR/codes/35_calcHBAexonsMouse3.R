#########################################################
## 35_calcHBAexonsMouse3.R
## Calculation of HBA along protein-coding genes (mouse)
## Rename unmasked mm10 genome: 20211215_nuCpos2_2.txt
## Download mm10 from NCBI (GRCm38_p6) and store in
## HOMEDIR/20211210_mouse_GRCm38_p6/NCBI_FASTA, and untar it.
library(Biostrings)
setwd(HOMEDIR)
setwd("20211210_mouse_GRCm38_p6/NCBI_FASTA")
setwd("genome_assemblies_genome_fasta/ncbi-genomes-2021-12-09")
mm10 <- readDNAStringSet(filepath = "GCF_000001635.26_GRCm38.p6_genomic.fna.gz")
setwd(HOMEDIR)
setwd("20211210_mouse_GRCm38_p6/NCBI_FASTA")
setwd("genome_assemblies_asm_struct")
setwd("ncbi-genomes-2021-12-10")
assemblyReport <- read.table(file = "GCF_000001635.26_GRCm38.p6_assembly_report.txt", skip = 43, header = FALSE)[, c(7, 10)]
mm10names <- names(mm10)
for(n in 1:length(mm10names)){
	NCBIname <- strsplit(mm10names[n], split = " ")[[1]][1]
	mm10names[n] <- assemblyReport$V10[assemblyReport$V7 == NCBIname]
}
names(mm10) <- mm10names
setwd(HOMEDIR)
setwd("20211210_mouse_GRCm38_p6")
setwd("RData")
save(mm10, file = "mm10.RData")


