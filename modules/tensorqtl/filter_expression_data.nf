params.run = true

process filter_expression_data {
    tag "${prefix}"
    queue 'normal'
    maxForks 2
    conda '/lustre/scratch118/humgen/resources/conda_envs/R3.6'
    memory = '55G'
    time '700m'
    cpus 1
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/filter_expression_data/", mode: 'symlink', pattern: "expression_data.filtered.bed", overwrite: true

    when:
    params.run

    input: 
    file(expression_data)
    val(variance_threshold)
    val(pcent_samples) 

    output: 
    tuple file("expression_data.filtered.bed"), emit: expression_data

    script:
    """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/R3.6/bin:\$PATH
exports R_LIBS=/lustre/scratch118/humgen/resources/rstudio_server_libs
cp ${expression_data} expression_data.filtered.bed
    """
}
