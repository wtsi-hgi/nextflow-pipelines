params.run = true

process 'baton_imeta' {
    tag "$sanger_sample_id"
    memory = '3G'
    time '120m'
    cpus 1
    maxForks 12
    errorStrategy { task.exitStatus in 5..6 ? 'ignore' : 'retry' }
    maxRetries 2
    publishDir "${params.outdir}/imeta/", mode: 'copy', pattern: "${sanger_sample_id}.imeta.tsv"

    when:
    params.run 

    input:
    val(sanger_sample_id)
    
    output:
    tuple val(sanger_sample_id), file("${sanger_sample_id}.imeta.tsv"), emit: imeta
    env(WORK_DIR), emit: work_dir_to_remove

    script:
    """
jq -n '{avus: [{attribute: \"sample\", value: \"${sanger_sample_id}\", o: \"=\"}, 
               {attribute: \"target\", value: \"1\", o: \"=\"}, 
               {attribute: \"manual_qc\", value: \"1\", o: \"=\"}]}' | \
singularity exec /software/hgi/containers/singularity-baton/baton.simg baton-metaquery \
  --zone seq \
  --obj --avu | \
jq --raw-output \
  '.[0].avus[] as \$t | [\$t.attribute, \$t.value] | @tsv' | \
sed s'/^sample\\t/sanger_sample_id\\t/' > to_filter.txt

cat to_filter.txt | \\
grep -P '^study\\t|^study_id\\t|^sanger_sample_id\\t|^study_id\\t|^manual_qc\\t|^id_run\\t|^lane\\t|^library_id\\t|^total_reads\\t|^library_type\\t' | \\
sort > to_transpose.txt

cat to_transpose.txt | cut -f 1 | tr '\\n' '\\t' > ${sanger_sample_id}.imeta.tsv
printf '\\n' >> ${sanger_sample_id}.imeta.tsv
cat to_transpose.txt | cut -f 2 | tr '\\n' '\\t' >> ${sanger_sample_id}.imeta.tsv

WORK_DIR=\$PWD
    """
}
// {attribute: \"manual_qc\", value: \"1\", o: \"=\"}, 

//singularity exec /software/hgi/containers/singularity-baton/baton.simg baton-metaquery \
