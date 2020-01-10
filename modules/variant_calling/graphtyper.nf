params.run = true

process graphtyper {
    memory '4G'
    tag "$bamlist_file"
    cpus 1
    disk '20 GB'
    scratch '/tmp'
    stageInMode 'copy'
    stageOutMode 'rsync'
    time '1000m'
    container "copy_number_v2"
    maxForks 60
    // containerOptions = "--bind /home/ubuntu"
    // errorStrategy 'terminate'
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    publishDir "${params.outdir}/graphtyper/", mode: 'symlink', overwrite: true, pattern: "commands_split.txt"
    maxRetries 2

    when:
    params.run
     
    input:
    file(bamlist_file)
    file(config_sh)
    file(bam_files)
    file(genome_fasta)   
    file(genome_fasta_fai)   

    output: 
    file("commands_split.txt"), emit: commands_split

    script:
""" 
cp -r /graphtyper-pipelines .
cd ./graphtyper-pipelines
bash make_graphtyper_pipeline.sh $bamlist_file $config_sh > ../commands_split.txt
"""
}
