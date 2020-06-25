params.run = true

process mosdepth_from_local {
    memory = '10G'
    time '240m'
    cpus 1
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    conda '/lustre/scratch118/humgen/resources/conda_envs/mosdepth'
    maxRetries 3
    maxForks 500
    publishDir "${params.outdir}/mosdepth/${study_alias}/${samplename}/", mode: 'copy', pattern: "${samplename}.mosdepth.*", overwrite: true

    when:
    params.run

    input: 
    set val(studyid), val(samplename), file(cram), val(study_alias), file(ref_fasta), file(capture_bed)

    output: 
    tuple val(studyid), val(samplename), file("${samplename}.mosdepth.region.dist.txt"), emit: study_samplename_regiondist
    tuple val(studyid), val(samplename), file("${samplename}.mosdepth.*"), emit: samplename_mosdepth
    //set val(samplename), file('*.cram') optional true // into ch_cram_files
    //file('*.lostcause.txt') optional true // into ch_lostcause_irods

    script:
    """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/mosdepth/bin:\$PATH
samtools index ${cram}

mosdepth --threads ${task.cpus} --fast-mode \\
--fasta ${ref_fasta} \\
--no-per-base \\
--by ${capture_bed} ${samplename} ${cram}

# clean-up
rm *.cram.crai
rm -f *.fa.fai
    """
}

// mosdepth --threads 4 --fast-mode --fasta ref.fasta --by noheaders.bed dd-africa8649298.out 33631_3#1.cram
// mosdepth --threads 4 --fast-mode --no-per-base --by 500 --fasta ref.fasta dd-africa8649298.out 33631_3#1.cram
