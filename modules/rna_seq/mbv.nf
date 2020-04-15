params.run = true 

process mbv {
    tag "mbv $samplename"
    memory = '10G'
    // conda '/lustre/scratch118/humgen/resources/conda/star'
    queue 'normal'
    time '700m'
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 2
    
    publishDir "${params.outdir}/mbv", mode: 'symlink'

    when:
    params.run

    input:
    set val(samplename), file(bam)
    file vcf_gz

    output:
    set val(samplename), file("${samplename}.bamstat.txt")
    set val(samplename), file("sorted_${samplename}_chr1.bam")

    script:
    """
samtools view -b $bam 1 > sorted_${samplename}_chr1.bam
samtools index sorted_${samplename}_chr1.bam

export PATH=/lustre/scratch118/humgen/resources/mbv/QTLtools_1.2_Ubuntu16.04_x86_64:\$PATH

QTLtools_1.2_Ubuntu16.04_x86_64 mbv --bam sorted_${samplename}_chr1.bam \\
--vcf $vcf_gz \\
--filter-mapping-quality 150 \\
--out ${samplename}.bamstat.txt
    """
}
