params.run = true

process prep_capture_bed {
    memory = '10G'
    time '240m'
    cpus 1
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    conda '/lustre/scratch118/humgen/resources/conda_envs/mosdepth'
    maxRetries 3
    maxForks 12
    publishDir "${params.outdir}/capture_bed/", mode: 'symlink', pattern: "${study_id}.capture.bed", overwrite: true

    when:
    params.run

    input: 
    set val(study_id), file(bed_file), val(remove_chr)

    output: 
    tuple val(study_id), file("${study_id}.capture.bed"), emit: study_capturebed

    script:
    """
if [ "$remove_chr" == "true" ]; then
    cat ${bed_file} | grep -v '^@' | grep -v '^chrY' | grep -v '^chrX' | grep -v '^X' | grep -v '^Y' | sed s/^chr//g > ${study_id}.capture.bed
else
    cat ${bed_file} | grep -v '^@' | grep -v '^chrY' | grep -v '^chrX' | grep -v '^X' | grep -v '^Y' | sed s'/^\\([1-9]\\)/chr\\1/'g > ${study_id}.capture.bed
fi
    """
}
