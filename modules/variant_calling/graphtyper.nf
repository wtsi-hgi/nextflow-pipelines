params.run = true

process graphtyper {
    memory '20G'
    tag "$graphtyper_command"
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
    tuple graphtyper_command, chr

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

$graphtyper_command
"""
}

