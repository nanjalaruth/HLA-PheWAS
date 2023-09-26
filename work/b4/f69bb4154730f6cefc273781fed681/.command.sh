#!/bin/bash -euo pipefail
Rscript /usr/local/bin/step1_fitNULLGLMM.R \
--plinkFile=merged_chr6_filtered_cauc_sex \
--phenoFile=cov_pheno.tsv \
--phenoCol=asthma_pheno \
--covarColList=sex,age,age_squared,age_sex \
--sampleIDColinphenoFile=IID \
--traitType=binary \
--outputPrefix=UKB_asthma_saige_out \
--nThreads=9 \
--LOCO=FALSE \
--IsOverwriteVarianceRatioFile=TRUE
