params.run = true

process vqsr_vcf {
    memory '75G'
    tag "vqsr $vcf"
    cpus 2
    //conda '/lustre/scratch118/humgen/hgi/projects/ibdx10/variant_calling/joint_calling/ibd_concat_nextflow/bcftools'
    //scratch '/tmp'
    //stageInMode 'copy'
    //stageOutMode 'copy'
    time '700m'
    queue 'normal'
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    publishDir "${params.outdir}/vqsr_vcf/$name/", mode: 'symlink', overwrite: true, pattern: "*.pdf"
    
    maxRetries 2

    when:
    params.run
     
    input:
    tuple file(name), file(vcf), file(tbi)
    
    output:
    tuple file("${name}.snps.tranches"), file("${name}.indels.tranches"), emit: tranches
    tuple file("${name}.snps.recal"), file("${name}.indels.recal"), emit: recal
    tuple file("*.R"), file("*.pdf"), emit: plots

    script:
    def simplename = vcf.getSimpleName()
""" 
sleep 10

CWD=\$PWD 
GATK=/lustre/scratch119/humgen/teams/hgi/users/ad7/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar
JAVA=/lustre/scratch119/realdata/mdt2/projects/interval_wgs/hgi/tools/java-8-openjdk-amd64/bin/java
ref_genome=/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/hs38DH.fa
dbSNP=/lustre/scratch118/humgen/resources/GATK/bundle/hg38/dbsnp_138.hg38.vcf.gz
hapmap=/lustre/scratch118/humgen/resources/GATK/bundle/hg38/hapmap_3.3.hg38.vcf.gz
omni_1000G=/lustre/scratch118/humgen/resources/GATK/bundle/hg38/1000G_omni2.5.hg38.vcf.gz
snps_1000G=/lustre/scratch118/humgen/resources/GATK/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
Mills_indel=/lustre/scratch118/humgen/resources/GATK/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
axiomPoly_resource=/lustre/scratch118/humgen/resources/GATK/bundle/hg38/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz

echo INDELS
singularity exec -B /lustre -B \$CWD -B /lustre/scratch118/humgen/resources /software/hgi/containers/gatk-4.1.0.0.simg /gatk/gatk --java-options "-XX:+UseSerialGC -Xmx64g -Xms64g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" VariantRecalibrator \
 	--reference \${ref_genome} \
 	-V ${vcf} \
 	--resource:mills,known=false,training=true,truth=true,prior=12.0 \${Mills_indel} \
 	--resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 \${axiomPoly_resource} \
 	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \${dbSNP} \
 	--trust-all-polymorphic \
 	-an QD \
 	-an MQRankSum \
 	-an ReadPosRankSum \
 	-an FS \
 	-an SOR \
     --use-allele-specific-annotations \
 	-mode INDEL \
 	--max-gaussians 4 \
 	-O  ${name}.indels.recal \
 	--tranches-file ${name}.indels.tranches \
 	--rscript-file ${name}.indels.plots.R \
 	-tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 \
 	-tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.1 \
 	-tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 \
 	-tranche 94.0 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 

echo SNPS
singularity exec -B /lustre -B \$CWD -B /lustre/scratch118/humgen/resources /software/hgi/containers/gatk-4.1.0.0.simg /gatk/gatk --java-options "-XX:+UseSerialGC -Xmx64g -Xms64g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" VariantRecalibrator \
	--reference \${ref_genome} \
	-V ${vcf} \
	--resource:hapmap,known=false,training=true,truth=true,prior=15.0 \${hapmap} \
	--resource:omni,known=false,training=true,truth=true,prior=12.0 \${omni_1000G} \
	--resource:1000G,known=false,training=true,truth=false,prior=10.0 \${snps_1000G} \
	--resource:dbsnp,known=true,training=false,truth=false,prior=7.0 \${dbSNP} \
	--trust-all-polymorphic \
	-an QD \
	-an MQ \
	-an MQRankSum \
	-an ReadPosRankSum \
	-an FS \
	-an SOR \
    --use-allele-specific-annotations \
	-mode SNP \
	--max-gaussians 6 \
	-O ${name}.snps.recal \
	--tranches-file ${name}.snps.tranches \
	--rscript-file ${name}.snps.plots.R \
	-tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 \
	-tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.1 \
	-tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 \
	-tranche 94.0 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 

"""
}

// /lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/output_vcf/stripped_vcf/VQSR_indel.sh ${vcf}
// /lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/output_vcf/stripped_vcf/VQSR_snp.sh ${vcf} && \\
// /lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/output_vcf/stripped_vcf/VQSR_indel_apply.sh ${vcf} && \\ 
// /lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/output_vcf/stripped_vcf/VQSR_snp_apply.sh ${vcf} \\
// bgzip ${simplename}.vqsr.vcf
// bcftools index ${simplename}.vqsr.vcf.gz
