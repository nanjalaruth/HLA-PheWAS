nextflow.enable.dsl=2

process saige_logistic_whole_1 {
    tag "Run logistic regression on ${dataset}_${subpop}"
    publishDir "${params.outdir}/Regression_results", mode: 'copy', overwrite: true
    label "bigmem"
    container "/apps/singularity/saige_1.1.6.3.sif"

    input:
        tuple val(dataset), path(bed), path(bim), path(fam), val(pheno_label), path(phenotype)

    output:
       tuple val(dataset), val(pheno_label), 
            file("${dataset}_${pheno_label}_saige_out.rda"),
               file("${dataset}_${pheno_label}_saige_out.varianceRatio.txt")

    script:
        out = "${dataset}_${pheno_label}_saige_out"
        base = bed.baseName

        """
        Rscript /usr/local/bin/step1_fitNULLGLMM.R \\
        --plinkFile=${base} \\
        --phenoFile=${phenotype} \\
        --phenoCol=${pheno_label}_pheno \\
        --covarColList=sex,baseline_age,Gansu,Haikou,Harbin,Henan,Hunan,Liuzhou,Qingdao,Sichuan,Suzhou,Zhejiang,hep_b_1_DPB1_Dosage,hep_b_2_DPB1_Dosage,hep_b_3_DPB1_Dosage,hep_b_4_DPB1_Dosage,hep_b_5_DPB1_Dosage,hep_b_1_DRB1_Dosage,hep_b_2_DRB1_Dosage,hep_b_3_DRB1_Dosage,hep_b_4_DRB1_Dosage,hep_b_5_DRB1_Dosage,hep_b_1_DQB1_Dosage,hep_b_2_DQB1_Dosage,hep_b_1_A_Dosage,hep_b_2_A_Dosage \\
        --sampleIDColinphenoFile=IID \\
        --traitType=binary \\
        --outputPrefix=${out} \\
        --nThreads=9 \\
        --LOCO=FALSE \\
        --IsOverwriteVarianceRatioFile=TRUE
        """
}

process saige_logistic_whole_2 {
    tag "Run logistic regression on ${dataset}_${subpop}"
    publishDir "${params.outdir}/Regression_results/SAIGE", mode: 'copy', overwrite: true
    label "bigmem"
    container "/apps/singularity/saige_1.1.6.3.sif"

    input:
        tuple val(dataset), path(vcf), path(index), path(id), val(pheno_label), path(rda), path(varianceRatio)

    output:
       tuple val(dataset), val(pheno_label), file(out)

    script:
        out = "${dataset}_${pheno_label}_binary.SAIGE.vcf.genotype.txt"

        """
        Rscript /usr/local/bin/step2_SPAtests.R \\
        --vcfFile=${vcf} \\
        --vcfFileIndex=${index} \\
        --vcfField=GT --chrom=6 \\
        --minMAF=0.01 --minMAC=1 \\
        --sampleFile=${id} \\
        --GMMATmodelFile=${rda} \\
        --varianceRatioFile=${varianceRatio} \\
        --SAIGEOutputFile=${out} \\
        --LOCO=FALSE

        """
}

process saige_logistic_1 {
    tag "Run logistic regression on ${dataset}_${subpop}"
    publishDir "${params.outdir}/Regression_results", mode: 'copy', overwrite: true
    label "bigmem"
    container "/apps/singularity/saige_1.1.6.3.sif"

    input:
        tuple val(dataset), path(bed), path(bim), path(fam), path(phenotype)

    output:
       tuple val(dataset), file("${dataset}_RA_saige_out.rda"), file("${dataset}_RA_saige_out.varianceRatio.txt")
    
    script:
        out = "${dataset}_RA_saige_out"
        base = bed.baseName

        """
        Rscript /usr/local/bin/step1_fitNULLGLMM.R \\
        --plinkFile=${base} --phenoFile=${phenotype} \\
        --phenoCol=RA_pheno \\
        --covarColList=sex,age,age_squared,age_sex \\
        --sampleIDColinphenoFile=IID --traitType=binary \\
        --outputPrefix=${out} --nThreads=9 \\
        --LOCO=FALSE --IsOverwriteVarianceRatioFile=TRUE

        """

}


process saige_logistic_2 {
    tag "Run logistic regression on ${dataset}_${subpop}"
    publishDir "${params.outdir}/Regression_results", mode: 'copy', overwrite: true
    label "bigmem"
    container "/apps/singularity/saige_1.1.6.3.sif"

    input:
        tuple val(dataset), path(vcf), path(index), path(id), path(rda), path(varianceRatio)

    output:
       tuple val(dataset), file(out)
    
    script:
        out = "${dataset}_RA_binary.SAIGE.vcf.genotype.txt"

        """
        Rscript /usr/local/bin/step2_SPAtests.R \\
        --vcfFile=${vcf} \\
        --vcfFileIndex=${index} \\
        --vcfField=GT --chrom=6 \\
        --minMAF=0.01 --minMAC=1 \\
        --sampleFile=${id} \\
        --GMMATmodelFile=${rda} \\
        --varianceRatioFile=${varianceRatio} \\
        --SAIGEOutputFile=${out} \\
        --LOCO=FALSE

        """
}

workflow{
     //step 0
     cov_pheno_ch = Channel.fromList(params.cov_pheno)
     saige_spop_ch = Channel.fromList(params.whole_ckb)
     input = saige_spop_ch.combine(cov_pheno_ch)
     saige_logistic_whole_1(input)
  
     //step0.1
    step1_out = saige_logistic_whole_1.out
    vcf_ch = Channel.fromList(params.whole_ckb_vcf)
    ids_ch = Channel.fromPath(params.whole_ids)
    reg_2_ch = vcf_ch.combine(ids_ch)
          .combine(step1_out, by:0)
    //reg_2_ch.view()
    saige_logistic_whole_2(reg_2_ch)

     //step 1
     geno_ch = Channel.fromList(params.geno)
     cov_ch = Channel.fromList(params.pheno)
     //Step 1.2 Run regression
     reg_sh = geno_ch     
	.combine(cov_ch, by:0)
    // reg_sh.view()
    //saige_logistic_1(reg_sh)


    //step2
    //step1_out = saige_logistic_1.out
    //vcf_ch = Channel.fromList(params.whole_ckb_vcf)
    //ids_ch = Channel.fromList(params.ids)
    //reg_2_ch = ids_ch.combine(vcf_ch)
      //    .combine(step1_out, by:0)
        //  .map{subpop, ids, dataset, vcf, index, rda, ratio 
	//	-> [subpop, vcf, index, ids, rda, ratio]}
    //reg_2_ch.view()
    //saige_logistic_2(reg_2_ch)

    //step3
    //manhattan_plot(saige_logistic_2.out)
}
