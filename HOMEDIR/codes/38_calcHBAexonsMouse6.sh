#########################################################
## 38_calcHBAexonsMouse6.sh
## Calculation of HBA along protein-coding genes (mouse)
## Preparation of unique nucleosomes
## liftover of nucleosomes: 20211214_nuCpos2_1.txt (2)	
## Store mm9ToMm10.over.chain.gz in 20211210_mouse_GRCm38_p6/liftOver_unique
conda activate
cd HOMEDIR
cd 20211210_mouse_GRCm38_p6/liftOver_unique
liftOver uniqueNuc_mm9.bed mm9ToMm10.over.chain.gz uniqueNuc_mm10.bed uniqueNuc_unlifted.bed


