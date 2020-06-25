params.run = true
params.remove_crams = true

process mosdepth_from_irods {
    memory = '10G'
    time '240m'
    cpus 1
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    conda '/lustre/scratch118/humgen/resources/conda_envs/mosdepth'
    maxRetries 3
    maxForks 12
    publishDir "${params.outdir}/mosdepth/${study_alias}/${samplename}/", mode: 'copy', pattern: "${samplename}.mosdepth.*", overwrite: true

    when:
    params.run

    input: 
    set val(studyid), val(samplename), val(study_alias), file(ref_fasta), file(capture_bed)

    output: 
    tuple val(studyid), val(samplename), file("${samplename}.mosdepth.region.dist.txt"), emit: study_samplename_regiondist
    tuple val(studyid), val(samplename), file("${samplename}.mosdepth.*"), emit: samplename_mosdepth
    //set val(samplename), file('*.cram') optional true // into ch_cram_files
    //file('*.lostcause.txt') optional true // into ch_lostcause_irods

    script:
    """
if bash -euo pipefail $workflow.projectDir/../bin/mosdepth/irods.sh -N ${task.cpus} -t ${studyid} -s ${samplename} ${params.dropqc}; then
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

# discard bam files, keep only cram
# rm -f *.bam

export PATH=/lustre/scratch118/humgen/resources/conda_envs/mosdepth/bin:\$PATH
# merge crams if more than one if found
n_crams=\$(find . -mindepth 1 -maxdepth 1 -type f -name "*.cram" -printf x | wc -c)
find . -mindepth 1 -maxdepth 1 -type f -name "*.cram"
if [ \$n_crams -eq 1 ];then
echo only 1 cram found
samtools index *.cram
mosdepth --threads ${task.cpus} --fast-mode \\
--fasta ${ref_fasta} \\
--no-per-base \\
--by ${capture_bed} ${samplename} *.cram
else
echo \$n_crams cram found
samtools merge \\
-@ ${task.cpus} -f merged.bam *.cram
samtools index merged.bam
mosdepth --threads ${task.cpus} --fast-mode \\
--fasta ${ref_fasta} \\
--no-per-base \\
--by ${capture_bed} ${samplename} merged.bam
fi

# clean-up
if [ "$params.remove_crams" == "true" ]; then
    echo removing bam and cram files
    rm -f *.cram
    rm -f *.cram.crai
    rm -f *.bam
    rm -f *.bam.bai
    rm -f *.fa.fai
fi
    """
}
