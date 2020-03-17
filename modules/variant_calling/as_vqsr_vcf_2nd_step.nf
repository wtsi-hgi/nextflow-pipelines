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
    tuple file(vcf), file(tbi)
    
    output:
    tuple file("${name}.snps.tranches"), file("${name}.indels.tranches"), emit: tranches
    tuple file("*.R"), file("*.pdf"), emit: plots

    script:
    def simplename = vcf.getSimpleName()

'''
echo indel_apply
singularity exec -B /lustre -B \$CWD -B /lustre/scratch118/humgen/resources -B /lustre/scratch118/humgen/hgi/projects/interval_wes/joint_calls/output_vcf/stripped_vcf /software/hgi/containers/gatk-4.1.0.0.simg /gatk/gatk --java-options "-XX:+UseSerialGC -Xmx64g -Xms64g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" ApplyVQSR \
 	-V ${vcf} \
 	-O recal_indel_${vcf} \
 	--tranches-file ${vcf}.indels.tranches \
    --recal-file ${vcf}.indels.recal \
 	--reference \${ref_genome} \
        --truth-sensitivity-filter-level 99.5 \
        -mode INDEL

echo SNP_apply
singularity exec -B /lustre -B \$CWD -B /lustre/scratch118/humgen/resources -B /lustre/scratch118/humgen/hgi/projects/interval_wes/joint_calls/output_vcf/stripped_vcf /software/hgi/containers/gatk-4.1.0.0.simg /gatk/gatk --java-options "-XX:+UseSerialGC -Xmx64g -Xms64g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" ApplyVQSR \
 	-V ${vcf} \
 	-O "recal_snp_${vcf} \
 	--tranches-file ${vcf}.snps.tranches  \
    --recal-file ${vcf}.snps.recal \
 	--reference \${ref_genome} \
        --truth-sensitivity-filter-level 99.7 \
        -mode SNP

'''

}