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
        if (filename ==~ /.*\.all_founds_in_irods\.txt/) "ils_logs/$filename"
        else if (filename ==~ /.*\.not_found\.txt/) "ils_missing/$filename"
    }
    publishDir "${params.outdir}/iget_cellranger/full_data/", mode: 'symlink', pattern: "cellranger_${samplename}"
    publishDir "${params.outdir}/iget_cellranger/raw_feature_bc_matrix/", mode: 'symlink', pattern: "cellranger_${samplename}/raw_feature_bc_matrix"
    publishDir "${params.outdir}/iget_cellranger/filtered_feature_bc_matrix/", mode: 'symlink', pattern: "cellranger_${samplename}/filtered_feature_bc_matrix"
    publishDir "${params.outdir}/iget_cellranger/metrics_summary/", mode: 'symlink', pattern: "cellranger_${samplename}/metrics_summary.csv"

    when:
    params.run 

    input:
    set val(samplename), val(run_id), val(sanger_sample_id)
    
    output:
    set val(samplename), file("cellranger_${samplename}") optional true
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


// iget -Kr /seq/31488/cellranger/cellranger302_count_31488_Pilot_study_of_dissociation_methods_for_human_gut_tissues8024873_GRCh38-3_0_0 .

//imeta qu -z seq -d study_id = ${study_id} and sample = ${sample} and target = 1 | grep collection | awk -F ' ' '{print \$2}' > collection.txt
//imeta qu -z seq -d study_id = ${study_id} and sample = ${sample} and target = 1 | grep dataObj | awk -F ' ' '{print \$2}' > dataObj.txt
//paste -d '/' collection.txt dataObj.txt > ${samplename}.${sample}.${study_id}.to_iget.txt

//sort -o ${samplename}.${sample}.${study_id}.to_iget.txt ${samplename}.${sample}.${study_id}.to_iget.txt
//num=1
//cat ${samplename}.${sample}.${study_id}.to_iget.txt | while read line
//do
//    if [ \$num -gt 1 ]
//    then
//        iget -K -f -v \${line} ${batch}.${samplename}.\${num}.cram
//        iget -K -f -v \${line}.crai ${batch}.${samplename}.\${num}.cram.crai
//    else
//        iget -K -f -v \${line} ${batch}.${samplename}.cram
//        iget -K -f -v \${line}.crai ${batch}.${samplename}.cram.crai
//    fi
//    ((num++))
//done
