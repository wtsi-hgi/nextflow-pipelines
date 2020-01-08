params.run = true

process vcf_split_all_chr {
    memory '4G'
    tag "$samplename"
    cpus 2
    disk '60 GB'
    scratch '/tmp'
    stageInMode 'copy'
    stageOutMode 'rsync'
    time '1000m'
    container "copy_number_v2"
    maxForks 40
    // containerOptions = "--bind /home/ubuntu"
    // errorStrategy 'terminate'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    publishDir "${params.outdir}/copy_number_splitchr/", mode: 'symlink', overwrite: true 
    maxRetries 3

    when:
    params.run
     
    input: 
    tuple val(samplename), file(cn_vcf)
    
    output: 
    tuple val(samplename), file("${samplename}.cn.chr_*.recode.vcf"), emit: samplename_cn_vcf

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
    for i in \{1..22\}
    do vcftools  --vcf  $cn_vcf  --chr \$i  --recode --recode-INFO-all --out  ${samplename}.cn.chr_\$i
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
