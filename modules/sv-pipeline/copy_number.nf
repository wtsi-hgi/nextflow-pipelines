params.run = true

process copy_number_v2 {
    memory '4G'
    tag "copy number $samplename"
    cpus 2
    disk '29 GB'
    scratch '/tmp'
    stageInMode 'copy'
    stageOutMode 'rsync'
    time '1000m'
    container "copY_number_v2"
    maxForks 50
    // containerOptions = "--bind /home/ubuntu"
    // errorStrategy 'terminate'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    publishDir "${params.outdir}/copy_number/", mode: 'symlink', overwrite: true 
    maxRetries 3

    when:
    params.run
     
    input: 
    tuple val(samplename), val(egan_id), file(hist_root_file), file(samplename_gt_vcf), file(samplename_cram_json)
    
    output: 
    tuple val(samplename), file("${samplename}.cn.vcf"), emit: samplename_cn_vcf
    file("${samplename}.cn.vcf"), emit: cn_vcf

    script:
    """ 
    export ROOTSYS=/root
    export PATH=/root/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/speedseq/bin:/miniconda/bin
    
    eval \"\$(conda shell.bash hook)\"
    conda activate py2

   create_coordinates \\
      -i ${samplename_gt_vcf} \\
      -o coordinates.txt

    svtools copynumber \\
      -i ${samplename_gt_vcf} \\
      -s ${samplename} \\
      --cnvnator cnvnator \\
      -w 100 \\
      -r ${hist_root_file} \\
      -c coordinates.txt \\
      > ${samplename}.cn.vcf
    """
}
