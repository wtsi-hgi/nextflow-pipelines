params.run = true

process vqsr_vcf_apply {
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
    tuple file(name), file(vcf), file(tbi), file(snp_recal), file(snp_recal_index), file(indel_recal), file(indel_recal_index), file(snp_tranch), file(indel_tranch)
    
    output:
    tuple file("recal_indel_${vcf}"), file("recal_snp_${vcf}"), emit: apply_output

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

echo indel_apply
singularity exec -B /lustre -B \$CWD -B /lustre/scratch118/humgen/resources -B /lustre/scratch118/humgen/hgi/projects/interval_wes/joint_calls/output_vcf/stripped_vcf /software/hgi/containers/gatk-4.1.0.0.simg /gatk/gatk --java-options "-XX:+UseSerialGC -Xmx64g -Xms64g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" ApplyVQSR \
 	-V ${vcf} \
 	-O recal_indel_${vcf} \
     -AS \
 	--tranches-file ${indel_tranch} \
    --recal-file ${indel_recal} \
 	--reference \${ref_genome} \
        --truth-sensitivity-filter-level 99.5 \
        -mode INDEL

echo SNP_apply
singularity exec -B /lustre -B \$CWD -B /lustre/scratch118/humgen/resources -B /lustre/scratch118/humgen/hgi/projects/interval_wes/joint_calls/output_vcf/stripped_vcf /software/hgi/containers/gatk-4.1.0.0.simg /gatk/gatk --java-options "-XX:+UseSerialGC -Xmx64g -Xms64g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" ApplyVQSR \
 	-V ${vcf} \
 	-O recal_snp_${vcf} \
     -AS \
 	--tranches-file ${snp_tranch}  \
    --recal-file ${snp_recal} \
 	--reference \${ref_genome} \
        --truth-sensitivity-filter-level 99.7 \
        -mode SNP

"""

}