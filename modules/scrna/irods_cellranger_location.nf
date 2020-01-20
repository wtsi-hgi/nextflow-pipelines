params.run = true

process 'iget_cellranger_location' {

    disk '100 GB'
    scratch '/tmp'
    stageInMode 'copy'
    stageOutMode 'rsync'
    
    cpus 1
    time '120m'
    
    tag "iget $samplename $location"
    memory = '3G'
    maxForks 12
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    publishDir "${params.outdir}/iget_cellranger/", mode: 'symlink',
        saveAs: { filename ->
        if (filename ==~ /.*\.all_founds_in_irods\.txt/) "ils_logs/${samplename}/$filename"
        else if (filename ==~ /.*\.not_found\.txt/) "ils_missing/${samplename}/$filename"
    }
    publishDir "${params.outdir}/iget_cellranger/full_data/", mode: 'symlink', pattern: "cellranger_${samplename}"
    publishDir "${params.outdir}/iget_cellranger/raw_feature_bc_matrix/", mode: 'symlink', pattern: "cellranger_${samplename}/raw_feature_bc_matrix"
    publishDir "${params.outdir}/iget_cellranger/filtered_feature_bc_matrix/", mode: 'symlink', pattern: "cellranger_${samplename}/filtered_feature_bc_matrix"
    publishDir "${params.outdir}/iget_cellranger/metrics_summary/", mode: 'symlink', pattern: "cellranger_${samplename}/metrics_summary.csv"
    publishDir "${params.outdir}/iget_cellranger/bams/", mode: 'symlink', pattern: "cellranger_${samplename}/*.bam"
    publishDir "${params.outdir}/iget_cellranger/bams/", mode: 'symlink', pattern: "cellranger_${samplename}/*.bai"
    publishDir "${params.outdir}/iget_cellranger/bams/", mode: 'symlink', pattern: "cellranger_${samplename}/raw_feature_bc_matrix/barcodes.tsv.gz"

    when:
    params.run 

    input:
    set val(samplename), val(location)
    
    output:
    tuple val(samplename), file("cellranger_${samplename}"), emit: cellranger_full_dir optional true
    tuple val(samplename), file("cellranger_${samplename}/*.bam"), file("cellranger_${samplename}/*.bam.bai"), file("cellranger_${samplename}/filtered_feature_bc_matrix/barcodes.tsv.gz"), emit: cellranger_sample_bam_barcodes optional true
    tuple val(samplename), file("cellranger_${samplename}/raw_feature_bc_matrix"), emit: cellranger_raw optional true
    tuple val(samplename), file("cellranger_${samplename}/filtered_feature_bc_matrix"), emit: cellranger_filtered optional true
    tuple val(samplename), file("cellranger_${samplename}/metrics_summary.csv"), emit: cellranger_metrics_summary optional true
    tuple val(samplename), file("$samplename}.all_founds_in_irods.txt"), file("${samplename}.not_found.txt"), emit: cellranger_missing optional true

  script:
    """
ils ${location} > ${samplename}.all_founds_in_irods.txt
ils ${location} | grep ${samplename} > found_in_irods.txt

iget -Kr ${location} cellranger_${samplename}
   """
}

