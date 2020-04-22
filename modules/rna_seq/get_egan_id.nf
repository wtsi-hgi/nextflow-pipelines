params.run = true 

process get_egan_id {
    tag "$samplename"
    memory = '4G'
    // conda '/lustre/scratch118/humgen/resources/conda/star'
    queue 'normal'
    time '100m'
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 2
    maxForks 50
    //publishDir "${params.outdir}/egan_id", mode: 'symlink'

    when:
    params.run

    input:
    set val(samplename), file(cram)

    output:
    tuple val(samplename), file(cram), file("${samplename}_egan_id.csv"), emit: samplename_egan_id_csv

    script:
    """
EGANID=\$(samtools view -H $cram | \\
    grep 'SM:EGAN' | \\
    sed s'/^.*SM:/SM:/'g | \\
    sed s'/SM://'g | \\
    awk '{print \$1;}')

echo samplename,egan_id > ${samplename}_egan_id.csv
echo ${samplename},\$EGANID >> ${samplename}_egan_id.csv
    """
}
