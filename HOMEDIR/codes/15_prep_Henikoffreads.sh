#########################################################
## 15_prep_Henikoffreads.sh
## Download Henikoff's reads GSE51949: 20170822_Fuse_Revise_2.txt
HOMEDIR=XXX
cd HOMEDIR
cd 20210721_GSE51949/SRA
cd SRX372005
mkdir fastq
for SNAME in SRR1022812 SRR1022815 SRR1022818 SRR1022821 SRR1022813 \
	SRR1022816 SRR1022819 SRR1022822 SRR1022814 SRR1022817 SRR1022820 SRR1022823
do
cd $SNAME
fastq-dump --split-files ${SNAME}.sra
mv ${SNAME}_1.fastq ../fastq
mv ${SNAME}_2.fastq ../fastq
cd ../
done
cd fastq
cat *_1.fastq > SRX372005_1.fastq
cat *_2.fastq > SRX372005_2.fastq
rm SRR*
cd HOMEDIR
cd 20210721_GSE51949/SRA
cd SRX372006
mkdir fastq
for SNAME in SRR1022824 SRR1022829 SRR1022834 SRR1022839 SRR1022844 \
	SRR1022825 SRR1022830 SRR1022835 SRR1022840 SRR1022845 \
	SRR1022826 SRR1022831 SRR1022836 SRR1022841 \
	SRR1022827 SRR1022832 SRR1022837 SRR1022842 \
	SRR1022828 SRR1022833 SRR1022838 SRR1022843
do
cd $SNAME
fastq-dump --split-files ${SNAME}.sra
mv ${SNAME}_1.fastq ../fastq
mv ${SNAME}_2.fastq ../fastq
cd ../
done
cd fastq
cat *_1.fastq > SRX372006_1.fastq
cat *_2.fastq > SRX372006_2.fastq
rm SRR*
cd HOMEDIR
cd 20210721_GSE51949/SRA
cd SRX372007
mkdir fastq
for SNAME in SRR1022846 SRR1022851 SRR1022856 SRR1022861 SRR1022866 \
	SRR1022847 SRR1022852 SRR1022857 SRR1022862 SRR1022867 \
	SRR1022848 SRR1022853 SRR1022858 SRR1022863 SRR1022868 \
	SRR1022849 SRR1022854 SRR1022859 SRR1022864 \
	SRR1022850 SRR1022855 SRR1022860 SRR1022865
do
cd $SNAME
fastq-dump --split-files ${SNAME}.sra
mv ${SNAME}_1.fastq ../fastq
mv ${SNAME}_2.fastq ../fastq
cd ../
done
cd fastq
cat *_1.fastq > SRX372007_1.fastq
cat *_2.fastq > SRX372007_2.fastq
rm SRR*
cd HOMEDIR
cd 20210721_GSE51949
mkdir bwa
for SNAME in SRX372005 SRX372006 SRX372007
do
bwa mem -t 8 reference_genome/sc_genome_2011.fasta \
SRA/$SNAME/fastq/${SNAME}_1.fastq SRA/$SNAME/fastq/${SNAME}_2.fastq > bwa/${SNAME}.sam
rm SRA/$SNAME/fastq/${SNAME}_1.fastq
rm SRA/$SNAME/fastq/${SNAME}_2.fastq
cd bwa
samtools view -u ${SNAME}.sam | samtools sort -o ${SNAME}.bam -
samtools index ${SNAME}.bam
samtools flagstat ${SNAME}.bam > ${SNAME}_flagstat.txt
samtools idxstats ${SNAME}.bam > ${SNAME}_idxstats.txt
cd ../
done


