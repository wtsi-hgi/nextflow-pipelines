params.run = true

process graphtyper_on_interval {
    memory '60G'
    tag "$chr $start $end"
    cpus 4
    // disk '20 GB'
    time '1400m'
    queue 'long'
    container "graphtyper"
    containerOptions = "--bind /lustre"
    // errorStrategy 'terminate'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    publishDir "${params.outdir}/graphtyper/", mode: 'symlink', overwrite: true, pattern: "graphtyper-pipelines/results/$chr/*.vcf.gz"
    publishDir "${params.outdir}/graphtyper/", mode: 'symlink', overwrite: true, pattern: "graphtyper-pipelines/results/$chr/*.vcf.gz.tbi"
    publishDir "${params.outdir}/graphtyper/", mode: 'symlink', overwrite: true, pattern: "graphtyper-pipelines/haps/$chr/*.vcf.gz"
    publishDir "${params.outdir}/graphtyper/", mode: 'symlink', overwrite: true, pattern: "graphtyper-pipelines/haps/$chr/*.vcf.gz.tbi"
    publishDir "${params.outdir}/graphtyper/", mode: 'symlink', overwrite: true, pattern: "graphtyper-pipelines/hap_calls/$chr/*.vcf.gz"
    publishDir "${params.outdir}/graphtyper/", mode: 'symlink', overwrite: true, pattern: "graphtyper-pipelines/hap_calls/$chr/*.vcf.gz.tbi"
    maxRetries 3

    when:
    params.run
     
    input:
    file(bamlist_file)
    file(config_sh)
    tuple chr, start, end

    output: 
    tuple file("graphtyper-pipelines/results/$chr/*.vcf.gz"),file("graphtyper-pipelines/results/$chr/*.vcf.gz.tbi"), emit: vcf
    tuple file("graphtyper-pipelines/haps/$chr/*.vcf.gz"),file("graphtyper-pipelines/haps/$chr/*.vcf.gz.tbi"), emit: haps_vcf
    tuple file("graphtyper-pipelines/hap_calls/$chr/*.vcf.gz"),file("graphtyper-pipelines/hap_calls/$chr/*.vcf.gz.tbi"), emit: hap_calls_vcf
    tuple stdout, emit: stdout

    script:
""" 
cp -r /graphtyper-pipelines .
cd ./graphtyper-pipelines
sed -i s'/ --no_sort//'g node_script.sh
rm -rf .git
mkdir ../tmp
cp call_script.sh ..

export CUSTOM_REGION_SIZE=\$(($end - $start + 1))
export CUSTOM_CHROMOSOMES=\"$chr\"
export CUSTOM_SLICE_SIZE=\$((CUSTOM_REGION_SIZE / 4))
export CUSTOM_PAD_SIZE=\$((CUSTOM_SLICE_SIZE / 4))

echo \$CUSTOM_REGION_SIZE
echo \$CUSTOM_CHROMOSOMES
echo \$CUSTOM_SLICE_SIZE
echo \$CUSTOM_PAD_SIZE

set -e; set -o pipefail; ./node_script.sh ../graphtyper_pipeline_config_on_interval.sh ../$bamlist_file $chr:$start
"""
}

