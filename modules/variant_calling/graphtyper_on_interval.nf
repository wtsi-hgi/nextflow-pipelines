params.run = true

process graphtyper_on_interval {
    tag "$chr $start $end"
    memory = {  10.GB + 10.GB * (task.attempt - 1) }
    cpus 4
    // disk '20 GB'
    scratch '/tmp'
    time '1400m'
    queue 'long'
    container "graphtyper"
    containerOptions = "--bind /lustre --bind /tmp"
    publishDir "${params.outdir}/graphtyper/", mode: 'symlink', overwrite: true, pattern: "./results/$chr/*.vcf.gz"
    publishDir "${params.outdir}/graphtyper/", mode: 'symlink', overwrite: true, pattern: "./results/$chr/*.vcf.gz.tbi"
    publishDir "${params.outdir}/graphtyper/", mode: 'symlink', overwrite: true, pattern: "./haps/$chr/*.vcf.gz"
    publishDir "${params.outdir}/graphtyper/", mode: 'symlink', overwrite: true, pattern: "./haps/$chr/*.vcf.gz.tbi"
    publishDir "${params.outdir}/graphtyper/", mode: 'symlink', overwrite: true, pattern: "./hap_calls/$chr/*.vcf.gz"
    publishDir "${params.outdir}/graphtyper/", mode: 'symlink', overwrite: true, pattern: "./hap_calls/$chr/*.vcf.gz.tbi"
    maxRetries 5
    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }

    when:
    params.run
     
    input:
    file(bamlist_file)
    file(config_sh_files)
    tuple chr, start, end

    output: 
    tuple file("./results/$chr/*.vcf.gz"),file("./results/$chr/*.vcf.gz.tbi"), emit: vcf
    tuple file("./haps/$chr/*.vcf.gz"),file("./haps/$chr/*.vcf.gz.tbi"), emit: haps_vcf
    tuple file("./hap_calls/$chr/*.vcf.gz"),file("./hap_calls/$chr/*.vcf.gz.tbi"), emit: hap_calls_vcf
    tuple stdout, emit: stdout

    script:
""" 
export CUSTOM_REGION_SIZE=\$(($end - $start + 1))
export CUSTOM_CHROMOSOMES=\"$chr\"
export CUSTOM_NUM_SLICES_RUNNING=1
export CUSTOM_SLICE_SIZE=\$((CUSTOM_REGION_SIZE))
export CUSTOM_PAD_SIZE=\$((CUSTOM_SLICE_SIZE / 4))

echo \$CUSTOM_REGION_SIZE
echo \$CUSTOM_CHROMOSOMES
echo \$CUSTOM_SLICE_SIZE
echo \$CUSTOM_PAD_SIZE
echo \$CUSTOM_NUM_SLICES_RUNNING

./node_script.sh ./graphtyper_pipeline_config_on_interval.sh ./$bamlist_file $chr:$start
"""
}

//cp -r /graphtyper-pipelines .
//cd ./graphtyper-pipelines
//sed -i s'/ --no_sort//'g node_script.sh
//rm -rf .git
//cp call_script.sh ..
