nextflow.preview.dsl=2
params.current_date = 'current_date'
params.run_solo = true

include solo from '../modules/solo/solo.nf' params(run: true, outdir: params.outdir)

workflow {

    Channel.fromPath("${baseDir}/../../inputs/filtered_h5.txt")
	.splitCsv(header: true)
	.map{row->tuple(row.samplename, file(row.filtered_h5))}
        .take(-1)
	.set{to_solo}
    //to_solo.view()

    if (params.run_solo) {
	solo(to_solo,
	      Channel.fromPath("${baseDir}/../../inputs/solo_params.json").collect())
    }

}
