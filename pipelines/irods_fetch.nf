nextflow.preview.dsl=2

// burkina faso study
// SEQCAP_WGS_GDAP_Burkino_Faso 2971

// ghana studies
// SEQCAP_WGS_GDAP_GHANA 3180
// GDAP_XTEN 3257
// GDAP_10x_Zulu_Ashanti 3779

params.run_from_studies = true // if true, will use study ids to find & fetch Irods cram files:
ch_input_studies = Channel.from('2971','3180','3257','3779')
// if false, will use samples listed in samples.tsv:
params.input_samples_csv = "${baseDir}/../../inputs/samples.tsv"
params.copy_mode = "move" // choose symlink, move or copy to stage in results dir

params.dropqc = ''
include baton_study from '../modules/irods_fetch/baton_study.nf' params(run: true, outdir: params.outdir, dropqc: params.dropqc)
include iget_sample_study from '../modules/irods_fetch/iget_sample_study.nf' params(run: true, outdir: params.outdir, copy_mode: params.copy_mode)
include iget_sample from '../modules/irods_fetch/iget_sample.nf' params(run: true, outdir: params.outdir, copy_mode: params.copy_mode)

workflow {
    
    // multiple study_ids allowed
    if (params.run_from_studies) {
	baton_study(ch_input_studies)
	
	to_iget = baton_study.out.samples_noduplicates_tsv
	    .map{a,b -> b}
	    .splitCsv(header: true, sep: '\t')
	    .map{row->tuple(row.sample, row.sample_supplier_name, row.study_id)}
	    .map{a,b,c-> tuple(a,c)}

	iget_sample_study(to_iget)
	
    }
    else {
	Channel.fromPath(params.input_samples_csv)
	    .splitCsv(header: true, sep: '\t')
	    .map { row -> row.sample }
	    .set{ch_to_iget}
	
	iget_sample(ch_to_iget)
    }
}
