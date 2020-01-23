params.run = true

process 'vireo' {
    tag "vireo $samplename"

    //// FCE
    disk '100 GB'
    scratch '/tmp'
    stageInMode 'symlink'
    stageOutMode 'rsync'
    cpus = 8
    time '8000m'
    tontainer "single_cell"
    containerOptions = "--bind /"
    memory = {  100.GB + 50.GB * (task.attempt-1) }
    ////// FCE 

    //// farm
    //cpus =   {  2 * Math.min(1, task.attempt) }
    //memory = {  30.GB + 20.GB * (task.attempt-1) }
    //time '300m'
    //container "single_cell"
    // queue 'basement'
    /// farm
    
    errorStrategy = { task.attempt <= 4 ? 'retry' : 'ignore' }
    maxRetries 4
    publishDir "${params.outdir}/vireo/", mode: 'symlink'

    when:
    params.run 

    input:
    set val(samplename), file(cell_data), val(n_pooled)
    
    output:
    tuple val(samplename), file("vireo_${samplename}"), emit: vireo_output_dir

  script:
   """
export PATH=/opt/conda/envs/conda_env/bin:/opt/conda/bin:\$PATH
vireo -c $cell_data -N $n_pooled -o vireo_${samplename}
   """
}
