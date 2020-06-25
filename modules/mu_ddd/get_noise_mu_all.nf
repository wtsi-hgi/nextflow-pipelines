params.run = true

process get_noise_mu_all {
    tag "$sample_sanger_id"
    memory = '8G'
    time '240m'
    cpus 1
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    conda '/lustre/scratch118/humgen/resources/conda_envs/R.4'
    maxRetries 1
    maxForks 12
    publishDir "${params.outdir}/noise_mu_all/", mode: 'copy', pattern: "${sample_sanger_id}.noise_mu_out.tsv", overwrite: true

    when:
    params.run

    input: 
    tuple val(sample_sanger_id)
    file(sample_for_noise_annotations)
    file(allDDD_crams_csv)
    file(noise_mu_rscript)

    output: 
    tuple val(sample_sanger_id), file("${sample_sanger_id}.noise_mu_out.tsv"), emit: noise_mu_out

    script:
    """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/R.4/bin:\$PATH
R_LIBS=/lustre/scratch118/humgen/resources/rlibs4.0.0

# get cram files
cat ${allDDD_crams_csv} |\\
  grep '${sample_sanger_id}' >> crams_to_get.txt

cat crams_to_get.txt | while read line || [[ -n \$line ]];
do
   iget -K -v \$line .
done

# index cram files:
find . -type f -name '*.cram' -exec samtools index {} \\;

# run R script:
Rscript ${noise_mu_rscript} \\
  ${sample_sanger_id} ${sample_for_noise_annotations} \\
  ${sample_sanger_id}.cram

# clean-up:
# rm ./*.cram
# rm ./*.cram.crai
# rm crams_to_get.txt
    """
}


//cat crams_to_get.txt | while read line || [[ -n \$line ]];
//do
//   cp \$line .
//done
