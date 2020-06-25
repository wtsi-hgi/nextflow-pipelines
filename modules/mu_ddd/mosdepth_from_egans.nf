params.run = true
params.remove_crams = true

process mosdepth_from_egans {
    memory = '10G'
    time '240m'
    cpus 1
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    conda '/lustre/scratch118/humgen/resources/conda_envs/mosdepth'
    maxRetries 3
    maxForks 12
    publishDir "${params.outdir}/mosdepth/${study_alias}/${samplename}/", mode: 'copy', pattern: "${samplename}.mosdepth.*", overwrite: true

    when:
    params.run

    input: 
    set val(studyid), val(samplename), val(study_alias), file(ref_fasta), file(capture_bed)

    output: 
    tuple val(studyid), val(samplename), file("${samplename}.mosdepth.region.dist.txt"), emit: study_samplename_regiondist
    tuple val(studyid), val(samplename), file("${samplename}.mosdepth.*"), emit: samplename_mosdepth
    //set val(samplename), file('*.cram') optional true // into ch_cram_files
    //file('*.lostcause.txt') optional true // into ch_lostcause_irods

    script:
    """
imeta qu -z seq -d sample_accession_number = \"${samplename}\" and target = \"library\" | grep collection | awk -F ' ' '{print \$2}' > collection.txt || true
imeta qu -z seq -d sample_accession_number = \"${samplename}\" and target = \"library\" | grep dataObj | awk -F ' ' '{print \$2}' > dataObj.txt || true
paste -d '/' collection.txt dataObj.txt | grep '.cram\$' > to_iget.txt || true
if [[ -f "to_iget.txt" && -s "to_iget.txt" ]]; then 
    echo "target library exist and not empty"
else 
    echo "target library not exist or empty" 
    echo try with target 1 instead
    imeta qu -z seq -d sample_accession_number = \"${samplename}\" and target = 1 | grep collection | awk -F ' ' '{print \$2}' > collection.txt || true
    imeta qu -z seq -d sample_accession_number = \"${samplename}\" and target = 1 | grep dataObj | awk -F ' ' '{print \$2}' > dataObj.txt || true
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
    iget -K -f -v \${line} \$filename
done


export PATH=/lustre/scratch118/humgen/resources/conda_envs/mosdepth/bin:\$PATH
# merge crams if more than one if found
n_crams=\$(find . -mindepth 1 -maxdepth 1 -type f -name "*.cram" -printf x | wc -c)
find . -mindepth 1 -maxdepth 1 -type f -name "*.cram"
if [ \$n_crams -eq 1 ];then
echo only 1 cram found
samtools index *.cram
mosdepth --threads ${task.cpus} --fast-mode \\
--fasta ${ref_fasta} \\
--no-per-base \\
--by ${capture_bed} ${samplename} *.cram
else
echo \$n_crams cram found
samtools merge \\
-@ ${task.cpus} -f merged.bam *.cram
samtools index merged.bam
mosdepth --threads ${task.cpus} --fast-mode \\
--fasta ${ref_fasta} \\
--no-per-base \\
--by ${capture_bed} ${samplename} merged.bam
fi

# clean-up
if [ "$params.remove_crams" == "true" ]; then
    echo removing bam and cram files
    rm -f *.cram
    rm -f *.cram.crai
    rm -f *.bam
    rm -f *.bam.bai
    rm -f *.fa.fai
    rm dataObj.txt 
    rm to_iget.txt
    rm collection.txt 
fi
    """
}
