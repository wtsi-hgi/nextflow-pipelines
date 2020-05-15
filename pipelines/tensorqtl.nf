nextflow.preview.dsl=2

params.run_tensorqtl = true

include tensorqtl from '../modules/tensorqtl/tensorqtl.nf' params(run:true, outdir: params.outdir)

workflow {

    if (params.run_tensorqtl) {

	// params.inputs_dir = "${baseDir}/../../inputs/test_data/"
	params.inputs_dir = "/lustre/scratch114/projects/ukbb_scrna/tensorqtl/test_data/"
	Channel.value(tuple(
	    "GEUVADIS.445_samples",
	    "${params.inputs_dir}GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered",
	    file("${params.inputs_dir}GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.bed"),
	    file("${params.inputs_dir}GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.bim"),
	    file("${params.inputs_dir}GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.fam"),
	    file("${params.inputs_dir}GEUVADIS.445_samples.expression.bed.gz"),
	    file("${params.inputs_dir}GEUVADIS.445_samples.covariates.txt")
	)).
	set{to_tensorqtl}

	to_tensorqtl.view()
	tensorqtl(to_tensorqtl)
    }
}
