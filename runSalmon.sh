#!/bin/bash
source ~/.bash_profile
for i in {1..36}
do
	sample="FY$i"
	sample2="/Volumes/IMAGES/RNASeq_Tfam_GFAP_Cre/$sample"
	echo ${sample}
	salmon quant -i ~/Dropbox/RaymondTools/cDNA/salmon/version15/mouse_release95 \
				 -l A \
				 -1 ${sample2}_R1.fastq.gz -2 ${sample2}_R2.fastq.gz \
				 -o salmonResult_V15/${sample}.quant \
				 --validateMappings \
				 --seqBias \
				 --useVBOpt \
				 --numBootstraps 30
done
