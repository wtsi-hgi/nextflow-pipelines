params.run = true

process 'iget' {
    tag "iget $sample"
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
    sample
    
  output:
    tuple sample, file("*.cram"), file ("*.crai"), emit: sample_cram_crai optional true

  script:
    """
imeta qu -z seq -d sample = ${sample} and target = 1 | grep collection | awk -F ' ' '{print \$2}' > collection.txt
imeta qu -z seq -d sample = ${sample} and target = 1 | grep dataObj | awk -F ' ' '{print \$2}' > dataObj.txt
paste -d '/' collection.txt dataObj.txt > ${sample}.to_iget.txt

sort -o ${sample}.to_iget.txt ${sample}.to_iget.txt
num=1
cat ${sample}.to_iget.txt | while read line
do
    if [ \$num -gt 1 ]
    then
        iget -K -f -v \${line} ${sample}.\${num}.cram
        iget -K -f -v \${line}.crai ${sample}.\${num}.cram.crai || true
    else
        iget -K -f -v \${line} ${sample}.cram
        iget -K -f -v \${line}.crai ${sample}.cram.crai || true
    fi
    ((num++))
done
   """
}
 // study_id = ${study_id} and
 // study_id = ${study_id} and
