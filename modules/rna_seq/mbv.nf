params.run = true 

process mbv {
    tag "mbv $samplename"
    memory = '10G'
    // conda '/lustre/scratch118/humgen/resources/conda/star'
    queue 'normal'
    time '700m'
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 2
    
    publishDir "${params.outdir}/mbv", mode: 'rellink'

    when:
    params.run

    input:
    set val(samplename), file(bam), file(bai)
    file vcf_gz
    file vcf_gz_csi

    output:
    set val(samplename), file("${samplename}.bamstat.txt")

    script:
    """
# requires samtools, bcftools, bgzip and qtltools
export PATH=/software/sciops/pkgg/bcftools/1.10.2/bin/bcftools:\$PATH
export PATH=/software/sciops/pkgg/samtools/1.10.0/bin/samtools:\$PATH
export PATH=/lustre/scratch118/humgen/resources/mbv/QTLtools_1.2_Ubuntu16.04_x86_64:\$PATH
export PATH=/software/hgi/installs/anaconda3/envs/nextflow20/bin/:\$PATH

# run qtltools mbv
QTLtools_1.2_Ubuntu16.04_x86_64 mbv --bam ${bam} \\
--vcf ${vcf_gz} \\
--filter-mapping-quality 150 \\
--out ${samplename}.bamstat.txt
    """
}
