params.run = true

params.ref_dir = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/"
params.bed_dir = "/lustre/scratch118/humgen/resources/exome/Homo_sapiens/"

process deepvariant {
    memory '8G'
    //tag "$cram_file"
    //cpus 1
    //disk '20 GB'
    //time '100m'
    //queue 'normal'
    container  = 'file:///software/hgi/containers/deepvariant_0.10_UKB.sif'
    containerOptions = "--bind /lustre --bind ${params.ref_dir}:/ref --bind ${params.bed_dir}:/bed_files --bind /tmp:/tmp"
    // errorStrategy 'terminate'
    //errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    publishDir "${params.outdir}/", mode: 'symlink', overwrite: true, pattern: "*gz*"
    maxRetries 3

    when:
    params.run
     
    input:
    tuple path(cram_file_sorted_dups_coord), path(cram_file_sorted_dups_coord_index)

    output:
    tuple path("${cram_file_sorted_dups_coord}.vcf.gz"), path("${cram_file_sorted_dups_coord}.vcf.gz.tbi"), path("${cram_file_sorted_dups_coord}.g.vcf.gz"), path("${cram_file_sorted_dups_coord}.g.vcf.gz.tbi")
    //tuple file("${cram_file}.sorted"), emit: indexes


    script:
""" 
/opt/deepvariant/bin/run_deepvariant --model_type=WES --customized_model=/opt/models/ukb_wes/model.ckpt-22236 --ref=/ref/hs38DH.fa --reads=${cram_file_sorted_dups_coord} --output_vcf=${cram_file_sorted_dups_coord}.vcf.gz --output_gvcf=${cram_file_sorted_dups_coord}.g.vcf.gz --intermediate_results_dir /tmp/deepvariant_tmp --num_shards=1 --regions=/bed_files/All_Exon_V5/GRCh38/S04380110_Padded_merged.bed
"""
}

