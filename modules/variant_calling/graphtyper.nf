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
    file(bamlist_file)
    file(config_sh)
    val(graphtyper_command)

    output: 
    tuple stdout, emit: stdout

    script:
""" 
cp -r /graphtyper-pipelines .
cd ./graphtyper-pipelines
rm -rf .git
mkdir ../tmp
cp call_script.sh ..

$graphtyper_command
"""
}

// mkdir ../tmp/graphtyper_calling.XXXXXX
