nextflow.preview.dsl = 2
// for_noise_annotations.nf

params.for_noise_annotations = "$baseDir/../../inputs/for_noise_annotations"
params.r_script_noise_mu = "$baseDir/../../inputs/get_noise_mu_all.R"

params.allDDD_irods_crams = "$baseDir/../../inputs/irods_cram_locations.txt"
// temporary, until cram files are actually moved from /lustre to Irods
//params.allDDD_lustre_crams = "$baseDir/../../inputs/allDDD_Source_CRAMs.csv"

params.run_get_noise_mu_all = true
include get_noise_mu_all from '../modules/mu_ddd/get_noise_mu_all.nf' params(run: true, outdir: params.outdir)

workflow {

    Channel.fromPath(params.for_noise_annotations)
	.splitCsv(header: false, sep: '\t')
	.map{row->tuple(row[0])}
	.unique()
        .set{samples_to_process}

    if (params.run_get_noise_mu_all) {
	get_noise_mu_all(samples_to_process.take(2),	
		    Channel.fromPath(params.for_noise_annotations).collect(),
		    Channel.fromPath(params.allDDD_irods_crams).collect(),
		    Channel.fromPath(params.r_script_noise_mu).collect())
    }
}

		    //Channel.fromPath(params.allDDD_lustre_crams).collect(),
