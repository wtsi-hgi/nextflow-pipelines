params.run = true

process solo {
    tag "${samplename}"
    queue 'gpu-normal'
    clusterOptions '-gpu num=1'
    // clusterOptions '-gpu "num=1:mode=exclusive_process"'
    maxForks 5
    container 'solo_and_scanpy'
    containerOptions = "--nv --bind /lustre/scratch118/humgen/resources/containers/solo_libs/:/home/ubuntu/ --bind /lustre --bind /tmp"
    memory = '4G'
    time '700m'
    cpus 1
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    publishDir "${params.outdir}/solo/", mode: 'symlink', pattern: "solo_$samplename", overwrite: true
    publishDir "${params.outdir}/solo_tsv/", mode: 'symlink', pattern: "${samplename}.tsv", overwrite: true
    publishDir "${params.outdir}/h5ad/$samplename/", mode: 'symlink', pattern: "${data_h5_file}ad", overwrite: true

    when:
    params.run

    input: 
    set val(samplename), file(data_h5_file)
    file(solo_params_json)

    output: 
    tuple val(samplename), file("solo_$samplename"), emit: sample_solooutdir
    tuple val(samplename), file("${samplename}.tsv"), emit: sample_solotsv
    tuple val(samplename), file("${data_h5_file}ad"), emit: sample_h5ad

    script:
    """
export PATH=/opt/conda/bin:/opt/conda/envs/solo/bin:\$PATH
eval \"\$(conda shell.bash hook)\"

conda activate scanpy
python $workflow.projectDir/../bin/solo/convert_h5_h5ad.py \\
-i $data_h5_file \\
-o ${data_h5_file}ad

conda deactivate 
conda activate solo
ls /home/ubuntu/
/home/ubuntu/solo $solo_params_json ${data_h5_file}ad -o solo_$samplename -g

conda deactivate 
conda activate scanpy
python $workflow.projectDir/../bin/solo/solo_outdir_to_tsv.py \\
-d $data_h5_file \\
-s solo_$samplename \\
-o ${samplename}.tsv

echo solo done
    """
}
// https://github.com/calico/solo
// clusterOptions "-gpu \"num=1:mode=exclusive_process\""
