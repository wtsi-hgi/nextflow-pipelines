params.run = true

process vep_vcf {
    memory '6G'
    tag "$name"
    cpus 1
    //conda '/lustre/scratch118/humgen/hgi/projects/ibdx10/variant_calling/joint_calling/ibd_concat_nextflow/bcftools'
    //scratch '/tmp'
    //stageInMode 'copy'
    //stageOutMode 'copy'
    time '700m'
    queue 'normal'
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    publishDir "${params.outdir}/vep_vcf/$name/", mode: 'symlink', overwrite: true, pattern: "*.vep.vcf.gz"
    publishDir "${params.outdir}/vep_vcf/$name/", mode: 'symlink', overwrite: true, pattern: "*.vep.vcf.gz.csi"
    
    maxRetries 2

    when:
    params.run
     
    input:
    tuple val(name), file(vcf), file(csi)
    
    output:
    tuple val(name), file("*.vep.vcf.gz"), file("*.vep.vcf.gz.csi"), emit: name_vcf_csi 

    script:
    def simplename = vcf.getSimpleName()
""" 
sleep 10
export DIR=\$PWD 
/software/singularity-v3.5.1/bin/singularity exec --bind /lustre --bind \$DIR --bind /software/hgi/containers/vep-loftee/mount_vep:/opt/vep/.vep --pwd /opt/vep/src/ensembl-vep /software/hgi/containers/vep-loftee/vep-loftee-light.img ./vep --verbose --offline --species homo_sapiens --assembly GRCh38 --everything --cache --dir_cache /opt/vep/.vep/ --dir_plugins /opt/vep/.vep/Plugins/ --species homo_sapiens --vcf --allele_number --force_overwrite --fasta /opt/vep/.vep/homo_sapiens/97_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --plugin LoF,loftee_path:/opt/vep/.vep/Plugins,human_ancestor_fa:/opt/vep/.vep/human_ancestor.fa.gz,gerp_bigwig:/opt/vep/.vep/gerp_conservation_scores.homo_sapiens.GRCh38.bw,gerp_database:/opt/vep/.vep/gerp_conservation_scores.homo_sapiens.GRCh38.bw,conservation_file:/opt/vep/.vep/loftee.sql,run_splice_predictions:0,donor_disruption_mes_cutoff:6,acceptor_disruption_mes_cutoff:7 -i \$DIR/$vcf -o \$DIR/${simplename}.vep.vcf

bgzip ${simplename}.vep.vcf
bcftools index ${simplename}.vep.vcf.gz
"""
}

