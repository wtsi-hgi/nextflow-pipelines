params.run = true

process rtg_vcfeval {
    tag "${sample}"
    memory = '4G'
    time '240m'
    cpus 1
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    maxForks 12
    publishDir "${params.outdir}/subset_sample/${sample}/", mode: 'symlink', overwrite: true
    conda '/lustre/scratch118/humgen/resources/conda_envs/rtg_tools'

    when:
    params.run

    input: 
    set file(vcf), file(tbi)
    val(sample)
    val(rtg_sample)
    set file(rtg_vcf), file(rtg_tbi)
    file(rtg_genome_template)
    file(rtg_region)
    val(tag)
    
    output: 
    tuple file("subset_sample.vcf.gz"), file('subset_sample.vcf.gz.tbi'), emit: vcf_tbi

    script:
    """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/rtg_tools/bin:\$PATH

rtg vcfeval \\
--baseline $rtg_vcf \\
--calls $vcf \\
--output ${tag}_${sample} \\
--evaluation-regions $bed \
--bed-regions $bed \
--template $template \
--sample "INTEGRATION,Sample_Diag-excap51-HG002-EEogPU" \
--threads=4 \
--output-mode "split"
    """
}
