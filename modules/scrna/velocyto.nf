params.run = true

process 'velocyto' {
    tag "velocyto $samplename"

    errorStrategy = 'ignore' // { task.attempt <= 4 ? 'retry' : 'ignore' }
    cpus 8 // {  2 * Math.min(1, task.attempt) }
    memory 90.GB // {  30.GB + 20.GB * (task.attempt-1) }
    maxRetries 0
    
    queue 'basement'
    time '30000m' // normal queue is 24 hours, long is 48h.
    
    //queue 'long'
    //time '2879m' // normal queue is 24 hours, long is 48h.
    maxForks 200
    
    conda "/lustre/scratch118/humgen/resources/conda_envs/velocyto"
    publishDir "${params.outdir}/velocyto/", mode: 'symlink'

    when:
    params.run 

    input:
    set val(samplename), file(cellranger_dir)
    file(cellranger_gtf)
    
    output:
    tuple val(samplename), file("velocyto_${samplename}"), emit: velocyto_output_dir

  script:
   """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/velocyto/bin:\$PATH
mkdir cellrangerdir && ln -s ../$cellranger_dir cellrangerdir/outs
velocyto run10x --verbose --samtools-threads 8 --samtools-memory 88000 cellrangerdir $cellranger_gtf
mv cellrangerdir/velocyto velocyto_${samplename}
rm -r cellrangerdir
   """
}
//eval \"\$(conda shell.bash hook)\"
//conda activate /lustre/scratch118/humgen/resources/conda_envs/velocyto
