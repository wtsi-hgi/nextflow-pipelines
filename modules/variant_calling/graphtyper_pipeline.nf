params.run = true

process graphtyper_pipeline {
    memory '4G'
    tag "$bamlist_file"
    cpus 1
    disk '20 GB'
    time '100m'
    container "graphtyper"
    maxForks 60
    containerOptions = "--bind /lustre"
    // errorStrategy 'terminate'
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    publishDir "${params.outdir}/graphtyper/", mode: 'symlink', overwrite: true, pattern: "commands_split.txt"
    maxRetries 2

    when:
    params.run
     
    input:
    file(bamlist_file)
    file(config_sh)

    output: 
    tuple file("commands_split.txt"), emit: commands_split

    script:
""" 
cp -r /graphtyper-pipelines .
cd ./graphtyper-pipelines
rm -rf .git
bash make_graphtyper_pipeline.sh ../$bamlist_file ../$config_sh > ../commands_split.txt
"""
}
