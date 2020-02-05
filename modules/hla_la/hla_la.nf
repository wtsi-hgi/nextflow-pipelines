params.run = true

process hla_la {
    memory '4G'
    tag "$samplename"
    cpus 22
    disk '80 GB'
    scratch '/tmp'
    stageInMode 'copy'
    stageOutMode 'rsync'
    time '1000m'
    container "hla-la-1.0.1"
    containerOptions = "--bind /tmp"
    maxForks 2
    // errorStrategy 'terminate'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    publishDir "${params.outdir}/hla_la/", mode: 'symlink', overwrite: true, pattern: "${samplename}.cn.chr1.recode.vcf"
    maxRetries 1

    when:
    params.run
     
    input: 
    tuple val(eganid), val(irods_cram)
    
    output: 
    tuple val(eganid), file("checksum.txt"), emit: out
    tuple val(eganid), file("${eganid}.cn.chr*.recode.vcf"), emit: out

    script:
    """ 
iget ${irods_cram} tmp.cram

CHKSUMIRODS=\$(ils -L $irods_cram | tail -n 1 | awk '{print \$1;}')
echo \$CHKSUMIRODS

echo \"(echo \$CHKSUMIRODS) tmp.cram\" | md5sum -c - > checksum.txt

s3cmd sync -r s3://hla/.hla-la/graphs/PRG_MHC_GRCh38_withIMGT /tmp/ # this is 29G could be put in /tmp?
HLA-LA.pl --BAM *.cram --graph /tmp/PRG_MHC_GRCh38_withIMGT --sampleID $eganid --workingDir . --maxThreads 22
    """
}
