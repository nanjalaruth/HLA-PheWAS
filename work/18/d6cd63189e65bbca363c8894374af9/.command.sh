#!/bin/bash -euo pipefail
Rscript /usr/local/bin/step2_SPAtests.R \
--vcfFile=merged_chr6_filtered_cauc_sex.vcf.gz \
--vcfFileIndex=merged_chr6_filtered_cauc_sex.vcf.gz.csi \
--vcfField=GT --chrom=6 \
--minMAF=0.01 --minMAC=1 \
--sampleFile=orig_ids.txt \
--GMMATmodelFile=UKB_asthma_saige_out.rda \
--varianceRatioFile=UKB_asthma_saige_out.varianceRatio.txt \
--SAIGEOutputFile=UKB_asthma_binary.SAIGE.vcf.genotype.txt \
--LOCO=FALSE
