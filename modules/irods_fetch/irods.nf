params.run = true

process 'iget' {
    tag "$sample"
    memory = '3G'
    time '120m'
    cpus 1
    maxForks 12
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    publishDir "${params.outdir}/iget/${sample}/", mode: 'symlink'

    when:
    params.run 

    input:
    val(sample)
    
  output:
    tuple sample, file("*.cram"), file ("*.crai"), emit: sample_cram_crai optional true

  script:
    """
imeta qu -z seq -d sample = ${sample} and target = \"library\" and manual_qc = 1 | grep collection | awk -F ' ' '{print \$2}' > collection.txt || true
imeta qu -z seq -d sample = ${sample} and target = \"library\" and manual_qc = 1 | grep dataObj | awk -F ' ' '{print \$2}' > dataObj.txt || true
paste -d '/' collection.txt dataObj.txt | grep '.cram\$' > to_iget.txt || true

if [[ -f "to_iget.txt" && -s "to_iget.txt" ]]; then 
    echo "target library exist and not empty"
else 
    echo "target library not exist or empty" 
    echo try with target 1 instead
    imeta qu -z seq -d sample = ${sample} and target = 1 and manual_qc = 1 | grep collection | awk -F ' ' '{print \$2}' > collection.txt || true
    imeta qu -z seq -d sample = ${sample} and target = 1 and manual_qc = 1 | grep dataObj | awk -F ' ' '{print \$2}' > dataObj.txt || true
    paste -d '/' collection.txt dataObj.txt | grep '.cram\$' > to_iget.txt || true
fi
if [[ -f "to_iget.txt" && -s "to_iget.txt" ]]; then 
    echo "target 1 exist and not empty"
else 
    echo "target 1 not exist or empty" 
    exit 1
fi
cat to_iget.txt | while read line
do
    filename=\$(echo \${line} | sed s'/^.*\\///'g)
    echo filename is \$filename
    iget --retries 2 -K -f -v \${line} \$filename
    iget --retries 2 -K -f -v \${line}.crai \$\{filename\}.crai
done
   """
}
