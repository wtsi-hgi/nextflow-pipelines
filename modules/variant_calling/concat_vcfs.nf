params.run = true

process concat_vcfs {
    memory '60G'
    tag "$vcfs_location $name"
    cpus 2
    // disk '20 GB'
    time '700m'
    queue 'normal'
    container "graphtyper"
    containerOptions = "--bind /lustre"
    // errorStrategy 'terminate'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    publishDir "${params.outdir}/vcfs_concat/", mode: 'symlink', overwrite: true, pattern: "${name}.vcf.gz"
    publishDir "${params.outdir}/vcfs_concat/", mode: 'symlink', overwrite: true, pattern: "${name}.vcf.gz.csi"
    maxRetries 3

    when:
    params.run
     
    input:
    val(vcfs_location)
    val(name)
    
    output: 
    tuple file("${name}.vcf.gz"), file("${name}.vcf.gz.csi"), emit: vcf_gz
    tuple file("to_concat.list"), file("to_concat_non_empty.list"), emit: concat_list

    script:
""" 
find $vcfs_location -name '*.vcf.gz' -exec sh -c \"echo {} \\\$(zcat {} | grep -v '^#' | wc -l)\" \\; >> to_concat.list

cat to_concat.list | grep -v 'gz 0\$' | sort > to_concat_non_empty.list

# bcftools concat -f to_concat_non_empty.list --allow-overlaps | bcftools sort -o ${name}.vcf.gz -O z
# bcftools index ${name}.vcf.gz
"""
}

