params.run = true

process graphtyper {
    memory '4G'
    tag "$graphtyper_command"
    cpus 4
    disk '20 GB'
    time '100m'
    container "graphtyper"
    maxForks 60
    containerOptions = "--bind /lustre"
    // errorStrategy 'terminate'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    publishDir "${params.outdir}/graphtyper_cmds/", mode: 'symlink', overwrite: true //, pattern: "commands_split.txt"
    maxRetries 3

    when:
    params.run
     
    input:
    val(graphtyper_command)

    output: 
    stdout, emit: stdout

    script:
""" 
cp -r /graphtyper-pipelines/* .
bash $graphtyper_command
"""
}
