## 01_prep_H3Q85C.sh
## Preparation of the H3-Q85C dataset: 20210706_nuCpos2_1.txt (1)
## Download a tar file containing bigwig (.bw) files 
## from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97290.
## Store the tar file in HOMEDIR/20210706_Chereji_GSE97290.
## Extract the downloaded tar file.
HOMEDIR=XXX
cd HOMEDIR
cd 20210706_Chereji_GSE97290
tar -xf GSE97290_RAW.tar

## Convert the .bw files into .wig files with ucsc-bigwigtowig (v377). 
bigWigToWig GSM2561057_Dyads_H3_CC_rep_1.bw GSM2561057_Dyads_H3_CC_rep_1.wig
bigWigToWig GSM2561058_Dyads_H3_CC_rep_2.bw GSM2561058_Dyads_H3_CC_rep_2.wig
bigWigToWig GSM2561059_Dyads_H3_CC_rep_3.bw GSM2561059_Dyads_H3_CC_rep_3.wig
