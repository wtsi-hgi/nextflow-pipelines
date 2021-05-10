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
    set val(samplename), file(bam)
    file vcf_gz

    output:
    set val(samplename), file("${samplename}.bamstat.txt")

    script:
    """
# requires samtools, bcftools, bgzip and qtltools
export PATH=/software/sciops/pkgg/bcftools/1.10.2/bin/bcftools:\$PATH
export PATH=/software/sciops/pkgg/samtools/1.10.0/bin/samtools:\$PATH
export PATH=/lustre/scratch118/humgen/resources/mbv/QTLtools_1.2_Ubuntu16.04_x86_64:\$PATH
export PATH=/software/hgi/installs/anaconda3/envs/nextflow20/bin/:\$PATH

# extract chromosome 1 from vcf
bcftools index $vcf_gz
bcftools view $vcf_gz --regions 1 > cohort.1.vcf.gz

# change chr naming convention to match cram files
echo "1 chr1" >> chr_name_conv.txt
bcftools annotate --rename-chrs chr_name_conv.txt cohort.1.vcf.gz \\
| bgzip > cohort.chr1.vcf.gz

# extract chromosome 1 from cram file
samtools index $bam
samtools view -b $bam chr1 > sorted_chr1.bam
samtools index sorted_chr1.bam

# run qtltools mbv
QTLtools_1.2_Ubuntu16.04_x86_64 mbv --bam sorted_chr1.bam \\
--vcf cohort.chr1.vcf.gz \\
--filter-mapping-quality 150 \\
--out ${samplename}.bamstat.txt

# clean-up intermediary files
rm ${bam}.crai 
rm ${vcf_gz}.csi
rm cohort.1.vcf.gz*
rm cohort.chr1.vcf.gz*
rm sorted_chr1.bam*
rm chr_name_conv.txt
    """
}
