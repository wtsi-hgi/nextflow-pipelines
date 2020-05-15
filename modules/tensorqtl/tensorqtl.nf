params.run = true

process tensorqtl {
    tag "${prefix}"
    queue 'gpu-normal'
    clusterOptions '-gpu num=1'
    // clusterOptions '-gpu "num=1:mode=exclusive_process"'
    maxForks 2
    conda '/lustre/scratch118/humgen/resources/conda_envs/singlecell_eqtl'
    memory = '55G'
    time '700m'
    cpus 1
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/tensorqtl/", mode: 'symlink', pattern: "${prefix}.cis_qtl.txt.gz", overwrite: true

    when:
    params.run

    input: 
    set val(prefix), val(plink_prefix_path), file(bed), file(bim), file(fam), file(expression_bed), file(covariates_file)

    output: 
    tuple val(prefix), file("${prefix}.cis_qtl.txt.gz"), emit: tensorqtl_parquet

    script:
    """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/singlecell_eqtl/bin:\$PATH

python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis
    """
}

// export PATH=/opt/conda/bin:/opt/conda/envs/tensorqtl/bin:\$PATH

// https://github.com/calico/tensorqtl
// clusterOptions "-gpu \"num=1:mode=exclusive_process\""

// container 'tensorqtl_and_scanpy'
// containerOptions = "--nv --bind /lustre/scratch118/humgen/resources/containers/tensorqtl_libs/:/home/ubuntu/ --bind /lustre --bind /tmp"
