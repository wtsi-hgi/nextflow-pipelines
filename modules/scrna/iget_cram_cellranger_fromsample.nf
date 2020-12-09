params.run = true

process 'iget_cram_cellranger' {
    tag "iget $sanger_sample_id"
    memory = '3G'
    time '120m'
    cpus 1
    maxForks 12
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    publishDir "${params.outdir}/iget_cellranger/", mode: 'symlink',
        saveAs: { filename ->
        if (filename ==~ /.*\.all_founds_in_irods\.txt/) "ils_logs/${sanger_sample_id}/$filename"
        else if (filename ==~ /.*\.not_found\.txt/) "ils_missing/${sanger_sample_id}/$filename"
        else if (filename ==~ /.*\.cellranger.\.not_found\.txt/) "cellranger_missing/${sanger_sample_id}/$filename"
    }
    publishDir "${params.outdir}/iget_cram_found/", mode: 'symlink', pattern: "*.cram"
    publishDir "${params.outdir}/iget_cram_found/", mode: 'symlink', pattern: "*.crai"
    publishDir "${params.outdir}/iget_cram_not_found/", mode: 'copy', pattern: "*.not_found.txt"
    
    publishDir "${params.outdir}/iget_cellranger/full_data/${biopsy_type}/${disease_status}/", mode: 'symlink', pattern: "${sanger_sample_id}"
    publishDir "${params.outdir}/iget_cellranger/raw_feature_bc_matrix/${biopsy_type}/${disease_status}/", mode: 'symlink', pattern: "${sanger_sample_id}/raw_feature_bc_matrix"
    publishDir "${params.outdir}/iget_cellranger/filtered_feature_bc_matrix/${biopsy_type}/${disease_status}/", mode: 'symlink', pattern: "${sanger_sample_id}/filtered_feature_bc_matrix"
    publishDir "${params.outdir}/iget_cellranger/metrics_summary//${biopsy_type}/${disease_status}/", mode: 'symlink', pattern: "${sanger_sample_id}/metrics_summary.csv"
    publishDir "${params.outdir}/iget_cellranger/bams/${biopsy_type}/${disease_status}/", mode: 'symlink', pattern: "${sanger_sample_id}/*.bam"
    publishDir "${params.outdir}/iget_cellranger/bams/${biopsy_type}/${disease_status}/", mode: 'symlink', pattern: "${sanger_sample_id}/*.bai"
    publishDir "${params.outdir}/iget_cellranger/bams/${biopsy_type}/${disease_status}/", mode: 'symlink', pattern: "${sanger_sample_id}/raw_feature_bc_matrix/barcodes.tsv.gz"

    when:
    params.run 

    input:
    set val(sanger_sample_id), val(biopsy_type), val(disease_status)
    
    output:
    env(WORK_DIR), emit: work_dir_to_remove
    tuple val(sanger_sample_id), file("*.cram"), emit: spname_batch_cram optional true
    tuple val(sanger_sample_id), file("${sanger_sample_id}.not_found.txt"), emit: iget_not_found optional true
    tuple val(sanger_sample_id), file("cellranger.${sanger_sample_id}.not_found.txt"), emit: cellranger_not_found optional true
    tuple val(sanger_sample_id), file("*.cram"), file ("*.crai"), emit: spname_batch_cram_crai optional true
    
    tuple val(sanger_sample_id), file("${sanger_sample_id}"), emit: cellranger_full_dir optional true
    tuple val(sanger_sample_id), file("${sanger_sample_id}/*.bam"), file("${sanger_sample_id}/*.bam.bai"), file("${sanger_sample_id}/raw_feature_bc_matrix/barcodes.tsv.gz"), emit: cellranger_sample_bam_barcodes optional true
    tuple val(sanger_sample_id), file("${sanger_sample_id}/raw_feature_bc_matrix"), emit: cellranger_raw optional true
    tuple val(sanger_sample_id), file("${sanger_sample_id}/filtered_feature_bc_matrix"), emit: cellranger_filtered optional true
    tuple val(sanger_sample_id), file("${sanger_sample_id}/metrics_summary.csv"), emit: cellranger_metrics_summary optional true
    tuple val(sanger_sample_id), file("${sanger_sample_id}.all_founds_in_irods.txt"), file("${sanger_sample_id}.not_found.txt"), emit: cellranger_missing optional true
    tuple val(sanger_sample_id), file("${sanger_sample_id}.sync_status.tsv"), emit: sync_status
    tuple val(sanger_sample_id), file("cellranger.iget.sh"), emit: iget_command optional true

    script:
    """
echo imeta collection
imeta qu -z seq -d sample = ${sanger_sample_id} and target = 1 and manual_qc = 1 | grep collection | awk -F ' ' '{print \$2}' > collection.txt || true
echo imeta collection done

if [ -s collection.txt ] 
then
    echo collection found
    imeta qu -z seq -d sample = ${sanger_sample_id} and target = 1 and manual_qc = 1 | grep dataObj | awk -F ' ' '{print \$2}' > dataObj.txt
    paste -d '/' collection.txt dataObj.txt >  ${sanger_sample_id}.to_iget.txt

    sort -o ${sanger_sample_id}.to_iget.txt ${sanger_sample_id}.to_iget.txt
    num=1
    cat ${sanger_sample_id}.to_iget.txt | while read line
    do
        echo num \$num
        if [ \$num -gt 1 ]
        then
            echo iget -K -f -v \${line} ${sanger_sample_id}.\${num}.cram
            echo iget -K -f -v \${line}.crai ${sanger_sample_id}.\${num}.cram.crai || true
        else
            echo iget -K -f -v \${line} ${sanger_sample_id}.cram
            echo iget -K -f -v \${line}.crai ${sanger_sample_id}.cram.crai || true
        fi
        ((num++))
    done

    export run_id=\$(head -n 1 collection.txt | sed s'/\\/seq\\///'g | sed s'/illumina\\/runs\\/..\\///'g | sed s'/\\/lane.*\$//'g)
    export run_id_2digits=\$(echo \${run_id} | head -c2)
    echo run_id is \$run_id
    echo run_id_2digits is \$run_id_2digits

    echo looking for cellranger data
    ils /seq/\${run_id}/cellranger/ > ${sanger_sample_id}.cellranger.all_founds_in_irods.txt || true
    ils /seq/\${run_id}/cellranger/ | grep ${sanger_sample_id} > cellranger.found_in_irods.txt || true
    ils /seq/illumina/\${run_id_2digits}/\${run_id}/cellranger/ >> ${sanger_sample_id}.cellranger.all_founds_in_irods.txt || true
    ils /seq/illumina/\${run_id_2digits}/\${run_id}/cellranger/ | grep ${sanger_sample_id} >> cellranger.found_in_irods.txt || true
    ils /seq/illumina/runs/\${run_id_2digits}/\${run_id}/cellranger/ >> ${sanger_sample_id}.cellranger.all_founds_in_irods.txt || true
    ils /seq/illumina/runs/\${run_id_2digits}/\${run_id}/cellranger/ | grep ${sanger_sample_id} >> cellranger.found_in_irods.txt || true
    ils /seq/illumina/cellranger/ >> ${sanger_sample_id}.cellranger.all_founds_in_irods.txt || true
    ils /seq/illumina/cellranger/ | grep ${sanger_sample_id} >> cellranger.found_in_irods.txt || true

    if [ -s cellranger.found_in_irods.txt ] 
    then 
            echo cellranger data found
            WORK_DIR=dont_remove
            cat cellranger.found_in_irods.txt  | awk '{print \$2}' | sed 's/^/iget -Kr /'g | sed 's/\$/ ${sanger_sample_id}/'g > cellranger.iget.sh
            bash cellranger.iget.sh
    else
            echo cellranger data not found
            WORK_DIR=\$PWD
            echo not found /seq/\${run_id}/cellranger/ grep ${sanger_sample_id} > ${sanger_sample_id}.cellranger.not_found.txt
    fi
    echo end cellranger fetch

else
    echo no collection found
    echo ${sanger_sample_id} > ${sanger_sample_id}.not_found.txt
    WORK_DIR=\$PWD
fi

if [ -f cellranger.iget.sh ]
then     
            printf 'sanger_sample_id\\tsync_status\\n${sanger_sample_id}\\tsynced' > ${sanger_sample_id}.sync_status.tsv
else
            printf 'sanger_sample_id\\tsync_status\\n${sanger_sample_id}\\tnot_synced' > ${sanger_sample_id}.sync_status.tsv
fi
echo done
   """
}

// export run_id=\$(head -n 1 collection.txt | sed 's/\\/seq\\///'g)
