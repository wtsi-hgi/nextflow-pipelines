params.run = true

process parent_trio {
    tag "$sample_sanger_id"
    memory = '8G'
    time '240m'
    cpus 1
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    conda '/lustre/scratch118/humgen/resources/conda_envs/R.4'
    maxRetries 1
    maxForks 12
    publishDir "${params.outdir}/parent_trio/", mode: 'copy', pattern: "${sample_sanger_id}.trio_bam2R.tsv", overwrite: true

    when:
    params.run

    input: 
    tuple val(sample_sanger_id), file(sample_file), val(father_sanger_id), val(mother_sanger_id)
    file(allDDD_crams_csv)
    file(parent_trio_rscript)

    output: 
    tuple val(sample_sanger_id), file("${sample_sanger_id}.trio_bam2R.tsv"), emit: sample_triobam2R

    script:
    """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/R.4/bin:\$PATH
R_LIBS=/lustre/scratch118/humgen/resources/rlibs4.0.0

# get trio cram files:
cat ${allDDD_crams_csv} |\\
  grep '${sample_sanger_id}\\|${father_sanger_id}\\|${mother_sanger_id}' >> crams_to_get.txt

cat crams_to_get.txt | while read line || [[ -n \$line ]];
do
   iget -K -v \$line .
done

# index trio cram files:
find . -type f -name '*.cram' -exec samtools index {} \\;

# run R script:
Rscript ${parent_trio_rscript} \\
  ${sample_sanger_id} ${sample_file} \\
  ${sample_sanger_id}.cram ${father_sanger_id}.cram ${mother_sanger_id}.cram 

# clean-up:
# rm ./*.cram
# rm ./*.cram.crai
# rm crams_to_get.txt
    """
}

