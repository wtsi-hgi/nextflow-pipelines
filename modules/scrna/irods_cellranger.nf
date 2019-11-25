params.run = true

process 'iget_cellranger' {
    tag "iget $samplename $run_id"
    memory = '3G'
    time '120m'
    cpus 1
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
    publishDir "${params.outdir}/iget_cellranger/bams/", mode: 'symlink', pattern: "cellranger_${samplename}/*.bam", saveAs: "${samplename}.bam"
    publishDir "${params.outdir}/iget_cellranger/bams/", mode: 'symlink', pattern: "cellranger_${samplename}/*.bai", saveAs: "${samplename}.bam.bai"
    publishDir "${params.outdir}/iget_cellranger/bams/", mode: 'symlink', pattern: "cellranger_${samplename}/raw_feature_bc_matrix/barcodes.tsv.gz"

    when:
    params.run 

    input:
    set val(samplename), val(run_id), val(sanger_sample_id)
    
    output:
    set val(samplename), file("cellranger_${samplename}") optional true
    set val(samplename), file("cellranger_${samplename}/*.bam"), file("cellranger_${samplename}/*.bam.bai"), file("cellranger_${samplename}/raw_feature_bc_matrix/barcodes.tsv.gz")  optional true
    set val(samplename), file("cellranger_${samplename}/raw_feature_bc_matrix") optional true
    set val(samplename), file("cellranger_${samplename}/filtered_feature_bc_matrix") optional true
    set val(samplename), file("cellranger_${samplename}/metrics_summary.csv") optional true
    set val(samplename), file("${sanger_sample_id}.all_founds_in_irods.txt"), file("${samplename}.not_found.txt") optional true

  script:
    """
ils /seq/${run_id}/cellranger/ > ${sanger_sample_id}.all_founds_in_irods.txt
ils /seq/${run_id}/cellranger/ | grep ${sanger_sample_id} > found_in_irods.txt

cat found_in_irods.txt  | awk '{print \$2}' | sed 's/^/iget -Kr /'g | sed 's/\$/ cellranger_${samplename}/'g > iget.sh
if [ -s found_in_irods.txt ] 
then
        bash iget.sh
else
        echo not found /seq/${run_id}/cellranger/ grep ${sanger_sample_id} > ${samplename}.not_found.txt
fi
   """
}

