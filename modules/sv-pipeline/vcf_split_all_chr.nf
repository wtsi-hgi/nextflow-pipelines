params.run = true

process vcf_split_all_chr {
    memory '4G'
    tag "$samplename"
    cpus 1
    disk '20 GB'
    scratch '/tmp'
    stageInMode 'copy'
    stageOutMode 'rsync'
    time '1000m'
    container "copy_number_v2"
    maxForks 40
    // containerOptions = "--bind /home/ubuntu"
    // errorStrategy 'terminate'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    publishDir "${params.outdir}/copy_number_splitchr/chr1/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr1.recode.vcf"
    publishDir "${params.outdir}/copy_number_splitchr/chr2/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr2.recode.vcf" 
    publishDir "${params.outdir}/copy_number_splitchr/chr3/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr3.recode.vcf" 
    publishDir "${params.outdir}/copy_number_splitchr/chr4/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr4.recode.vcf" 
    publishDir "${params.outdir}/copy_number_splitchr/chr5/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr5.recode.vcf" 
    publishDir "${params.outdir}/copy_number_splitchr/chr6/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr6.recode.vcf" 
    publishDir "${params.outdir}/copy_number_splitchr/chr7/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr7.recode.vcf" 
    publishDir "${params.outdir}/copy_number_splitchr/chr8/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr8.recode.vcf" 
    publishDir "${params.outdir}/copy_number_splitchr/chr9/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr9.recode.vcf" 
    publishDir "${params.outdir}/copy_number_splitchr/chr10/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr10.recode.vcf"
    publishDir "${params.outdir}/copy_number_splitchr/chr11/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr11.recode.vcf" 
    publishDir "${params.outdir}/copy_number_splitchr/chr12/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr12.recode.vcf"
    publishDir "${params.outdir}/copy_number_splitchr/chr13/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr13.recode.vcf" 
    publishDir "${params.outdir}/copy_number_splitchr/chr14/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr14.recode.vcf" 
    publishDir "${params.outdir}/copy_number_splitchr/chr15/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr15.recode.vcf" 
    publishDir "${params.outdir}/copy_number_splitchr/chr16/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr16.recode.vcf" 
    publishDir "${params.outdir}/copy_number_splitchr/chr17/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr17.recode.vcf" 
    publishDir "${params.outdir}/copy_number_splitchr/chr18/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr18.recode.vcf" 
    publishDir "${params.outdir}/copy_number_splitchr/chr19/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr19.recode.vcf" 
    publishDir "${params.outdir}/copy_number_splitchr/chr20/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr20.recode.vcf" 
    publishDir "${params.outdir}/copy_number_splitchr/chr21/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr21.recode.vcf" 
    publishDir "${params.outdir}/copy_number_splitchr/chr22/", mode: 'symlink', overwrite: true, pattern "${params.outdir}/copy_number_splitchr/chr1/${samplename}.cn.chr22.recode.vcf" 
    maxRetries 3

    when:
    params.run
     
    input: 
    tuple val(samplename), file(cn_vcf)
    
    output: 
    tuple val(samplename), file("${samplename}.cn.chr*.recode.vcf"), emit: samplename_cn_vcf

    script:
    """ 
export ROOTSYS=/root
export MANPATH=/root/man:/usr/local/man:/usr/local/share/man:/usr/share/man
export USER_PATH=/home/ubuntu/error/speedseq/bin/:/home/ubuntu/anaconda3/envs/py2/bin:/home/ubuntu/anaconda3/condabin:/usr/local/go/bin:/home/ubuntu/error/root/bin:/usr/local/go/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/home/ubuntu/go/bin:/home/ubuntu/go/bin:/bin:/usr/bin:/sbin:/usr/sbin:/usr/local/bin:/usr/local/sbin
export LD_LIBRARY_PATH=/root/lib:/.singularity.d/libs
export LIBPATH=/root/lib
export JUPYTER_PATH=/root/etc/notebook
export DYLD_LIBRARY_PATH=/root/lib
export PYTHONPATH=/root/lib
export SHLIB_PATH=/root/lib
export CMAKE_PREFIX_PATH=/root
export CLING_STANDARD_PCH=none

export PATH=/root/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/speedseq/bin:/miniconda/bin:\$PATH

    eval \"\$(conda shell.bash hook)\"
    conda activate py2
    for i in {1..22}
        do vcftools  --vcf  $cn_vcf  --chr chr\$i  --recode --recode-INFO-all --out  ${samplename}.cn.chr\$i
    done
    """
}

//vcftools --vcf ${samplename_gt_vcf} --not-chr chrX --not-chr chrY --recode --recode-INFO-all --out ${samplename}.noXY

//   create_coordinates \\
//      -i ${samplename_gt_vcf} \\
//      -o coordinates.txt
// cat coordinates.txt | grep -v chrX | grep -v chrY > coordinates_nochrXY.txt
   
//   svtools copynumber \\
//      -i ${samplename}.noXY.recode.vcf \\
//      -s ${egan_id} \\
//      --cnvnator cnvnator \\
//      -w 100 \\
//      -r ${hist_root_file} \\
//      -c coordinates_nochrXY.txt \\
//      > ${samplename}.noXY.cn.vcf

      // -i ${samplename_gt_vcf} \\
