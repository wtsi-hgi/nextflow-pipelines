params.run = true

process 'iget_crams' {
    tag "iget $samplename $study_id"
    memory = '3G'
    time '120m'
    cpus 1
    maxForks 12
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    publishDir "${params.outdir}/iget_found/", mode: 'symlink', pattern: "*.cram"
    publishDir "${params.outdir}/iget_found/", mode: 'symlink', pattern: "*.crai"
    publishDir "${params.outdir}/iget_not_found/", mode: 'copy', pattern: "*.not_found.txt"

    when:
    params.run 

    input:
    set val(samplename), val(batch), val(sample), val(study_id)
    
  output:
    tuple val(samplename), val(batch), file("*.cram"), emit: spname_batch_cram optional true
    tuple val(samplename), file("${samplename}.${batch}.${sample}.${study_id}.not_found.txt"), emit: iget_not_found optional true
    tuple val(samplename), val(batch), file("*.cram"), file ("*.crai"), emit: spname_batch_cram_crai optional true

  script:
    """
echo imeta collection
imeta qu -z seq -d study_id = ${study_id} and sample = ${sample} and target = 1 | grep collection | awk -F ' ' '{print \$2}' > collection.txt || true
echo imeta collection done

if [ -s collection.txt ] 
then
    echo collection found
    imeta qu -z seq -d study_id = ${study_id} and sample = ${sample} and target = 1 | grep dataObj | awk -F ' ' '{print \$2}' > dataObj.txt
    paste -d '/' collection.txt dataObj.txt > ${samplename}.${sample}.${study_id}.to_iget.txt

    sort -o ${samplename}.${sample}.${study_id}.to_iget.txt ${samplename}.${sample}.${study_id}.to_iget.txt
    num=1
    cat ${samplename}.${sample}.${study_id}.to_iget.txt | while read line
    do
        echo num \$num
        if [ \$num -gt 1 ]
        then
            iget -K -f -v \${line} ${batch}.${samplename}.\${num}.cram
            iget -K -f -v \${line}.crai ${batch}.${samplename}.\${num}.cram.crai || true
        else
            iget -K -f -v \${line} ${batch}.${samplename}.cram
            iget -K -f -v \${line}.crai ${batch}.${samplename}.cram.crai || true
        fi
        ((num++))
    done
else
    echo no collection found
    echo ${samplename},${batch},${sample},${study_id} > ${samplename}.${batch}.${sample}.${study_id}.not_found.txt
fi
echo done
   """
}
