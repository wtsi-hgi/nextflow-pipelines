params.run = true

process sect_concat_vcfs {
    memory '6G'
    tag "$batch"
    cpus 1
    conda '/lustre/scratch118/humgen/hgi/projects/ibdx10/variant_calling/joint_calling/ibd_concat_nextflow/bcftools'
    scratch '/tmp'
    stageInMode 'copy'
    stageOutMode 'copy'
    time '700m'
    queue 'normal'
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    publishDir "${params.outdir}/concat/$batch/", mode: 'symlink', overwrite: true, pattern: "*.sorted.vcf.gz"
    publishDir "${params.outdir}/concat/$batch/", mode: 'symlink', overwrite: true, pattern: "*.sorted.vcf.gz.csi"
    publishDir "${params.outdir}/concat/$batch/", mode: 'symlink', overwrite: true, pattern: "vcf_files_sorted"
    maxRetries 2

    when:
    params.run
     
    input:
    tuple val(batch), file(vcf_files)
    file(bed)
    
    output:
    tuple val(batch), file("*.sorted.vcf.gz"), file("*.sorted.vcf.gz.csi"), emit: batch_vcf 
    tuple val(batch), file("vcf_files_sorted"), emit: batch_vcflist 

    script:
""" 
for vcf_file in *.output.vcf.gz
do
    bcftools norm -O z --rm-dup all -R ${bed} -o intersected.\${vcf_file} \${vcf_file}
    bcftools index intersected.\${vcf_file}
done

ls *.output.vcf.gz | grep -v tbi | grep -v csi | grep intersected > vcf_files
cat vcf_files | sort > vcf_files_sorted

bcftools concat -n -O z -f vcf_files_sorted -o batch.vcf.gz
export FIRST=\$(cat vcf_files_sorted | head -n 1 | sed s'/intersected.//'g | sed s'/.output.vcf.gz//'g)
bcftools sort -O z batch.vcf.gz -o \$FIRST.sorted.vcf.gz
bcftools index \$FIRST.sorted.vcf.gz
"""
}

