params.run = true

process 'iget_crams' {
    tag "iget $samplename"
    memory = '3G'
    time '120m'
    cpus 1
    maxForks 12
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    publishDir "${params.outdir}/iget/", mode: 'symlink'

    when:
    params.run 

    input:
    set val(samplename), val(batch), val(sample), val(study_id)
    
  output:
    set val(samplename), val(batch), file("*.cram"), file ("*.crai")

  script:
    """
imeta qu -z seq -d study_id = ${study_id} and sample = ${sample} and target = 1 | grep collection | awk -F ' ' '{print \$2}' > collection.txt
imeta qu -z seq -d study_id = ${study_id} and sample = ${sample} and target = 1 | grep dataObj | awk -F ' ' '{print \$2}' > dataObj.txt
paste -d '/' collection.txt dataObj.txt > ${samplename}.${sample}.${study_id}.to_iget.txt

sort -o ${samplename}.${sample}.${study_id}.to_iget.txt ${samplename}.${sample}.${study_id}.to_iget.txt
num=1
cat ${samplename}.${sample}.${study_id}.to_iget.txt | while read line
do
    if [ \$num -gt 1 ]
    then
        iget -K -f -v \${line} ${batch}.${samplename}.\${num}.cram
        iget -K -f -v \${line}.crai ${batch}.${samplename}.\${num}.cram.crai
    else
        iget -K -f -v \${line} ${batch}.${samplename}.cram
        iget -K -f -v \${line}.crai ${batch}.${samplename}.cram.crai
    fi
    ((num++))
done
   """
}
