#########################################################
## 33_calcHBAexonsMouse1.sh
## Calculation of HBA along protein-coding genes (mouse)
## Liftover of nucleosomes: 20170912_mouse_liftOver_1.txt (1)
## Download and untar GSE82127_RAW.tar
## from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82127
## Download chain.gz and store it in GSE82127_RAW
## from http://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver
cd HOMEDIR
cd GSE82127_RAW
for CHRNUM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y
do
CHR=chr${CHRNUM}
awk '{{START = $2 - 1}{print ($1, START, $2, $3)}}' \
${CHR}_Chemical_NCPscore.txt > ${CHR}_NCPscore_mm9.bed
liftOver ${CHR}_NCPscore_mm9.bed mm9ToMm10.over.chain.gz ${CHR}_NCPscore_mm10.bed unlifted_${CHR}_NCPscore_mm10.bed
awk 'BEGIN{FS=" ";OFS="    "}{print ($1, $3, $4)}' ${CHR}_NCPscore_mm10.bed > ${CHR}_Chemical_NCPscore_mm10.txt
rm ${CHR}_NCPscore_mm9.bed
rm ${CHR}_NCPscore_mm10.bed
awk 'BEGIN{FS=" ";OFS="    "}{if($4 > 0){print ($1, $3, $4)}}' unlifted_${CHR}_NCPscore_mm10.bed \
> unlifted_${CHR}_NCPscore_mm10.txt
rm unlifted_${CHR}_NCPscore_mm10.bed
done


