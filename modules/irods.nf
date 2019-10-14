params.run = true

process 'iget' {
    tag "iget $samplename"
    memory = '3G'
    time '120m'
    cpus 1
    maxForks 12
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    publishDir "${params.outdir}/iget/", mode: 'symlink'
    publishDir "${params.outdir}/iget_move/", mode: 'move'

    when:
    params.run 

    input:
    set val(samplename), val(sample), val(study_id)
    
  output:
    set val(samplename), file("*.cram"), file ("*.crai") optional true

  script:
    """
imeta qu -z seq -d study_id = ${study_id} and sample = ${sample} and target = 1 | grep collection | awk -F ' ' '{print \$2}' > collection.txt
imeta qu -z seq -d study_id = ${study_id} and sample = ${sample} and target = 1 | grep dataObj | awk -F ' ' '{print \$2}' > dataObj.txt
paste -d '/' collection.txt dataObj.txt > ${samplename}.${sample}.${study_id}.to_iget.txt

cat ${samplename}.${sample}.${study_id}.to_iget.txt | while read line
do
   iget -K -f -v \${line} ${samplename}.cram
   iget -K -f -v \${line}.crai ${samplename}.cram.crai
done
   """
}
