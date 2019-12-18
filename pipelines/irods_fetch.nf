nextflow.preview.dsl=2

include iget from '../modules/irods_fetch/irods.nf' params(run: true, outdir: params.outdir)

workflow {

    Channel.fromPath("${baseDir}/../../inputs/samples.tsv")
	.splitCsv(header: true, sep: '\t')
	.map { row -> tuple(row.sample) }
	.set{ch_to_iget}

    iget(ch_to_iget)
}
