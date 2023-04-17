#########################################################
## 16_prep_Henikoffdata1.sh
## Preparation of Henikoff's data: 20210721_nuCpos2_8.txt (1)
HOMEDIR=XXX
cd HOMEDIR
cd 20210721_GSE51949
mkdir bedGraph
for SNAME in SRX372005 SRX372006 SRX372007
do
bedtools genomecov -ibam bwa/${SNAME}.bam \
-d -5 -strand + > bedGraph/${SNAME}_Watson.bedGraph
bedtools genomecov -ibam bwa/${SNAME}.bam \
-d -5 -strand - > bedGraph/${SNAME}_Crick.bedGraph
done


