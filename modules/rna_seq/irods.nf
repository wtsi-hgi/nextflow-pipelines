params.run = true
params.dropqc = ""

process iget_cram {
    tag "iget cram ${samplename} ${studyid}"
    memory = '10G'
    time '240m'
    cpus 1
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }
    maxRetries 4
    maxForks 12
    publishDir "${params.outdir}/irods_lost/${samplename}/", mode: 'symlink', pattern: "*.lostcause.txt", overwrite: true
    publishDir "${params.outdir}/irods_crams/${samplename}/", mode: 'symlink', pattern: "*.cram", overwrite: true
    publishDir "${params.outdir}/irods_crams/${samplename}/", mode: 'symlink', pattern: "*.crai", overwrite: true

    when:
    params.run

    input: 
    val samplename //from sample_list_irods.flatMap{ it.readLines() }
    val studyid //from sample_list_irods.flatMap{ it.readLines() }

    output: 
    set val(samplename), file('*.cram') optional true // into ch_cram_files
    file('*.lostcause.txt') optional true // into ch_lostcause_irods
    set val(samplename), file('*.crai') optional true // into ch_cram_files

    script:
    """
    if bash -euo pipefail $workflow.projectDir/../bin/rna_seq/irods.sh -N ${task.cpus} -t ${studyid} -s ${samplename} ${params.dropqc}; then
      true
    else
      stat=\$?
      if [[ \$stat == 64 ]];
        then tag='nofiles';
        echo -e "${samplename}\\tirods\\t\$tag" > ${samplename}.lostcause.txt
      else          
        tag='UNKNOWN'
        echo -e "${samplename}\\tirods\\t\$tag" > ${samplename}.lostcause.txt
        exit \$stat
      fi
    fi
    """
}
