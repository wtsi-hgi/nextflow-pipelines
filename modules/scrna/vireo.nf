params.run = true

process 'vireo' {
    tag "vireo $samplename"
    container "single_cell"
    memory = '3G'
    time '120m'
    cpus 1
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/vireo/", mode: 'symlink'

    when:
    params.run 

    input:
    set val(samplename), file(cell_data), val(n_donor)
    
    output:
    set val(samplename), file("cellsnp_${samplename}") 

  script:
   """
vireo -c $cell_data -N $n_donor -o vireo_${samplename}
   """
}
