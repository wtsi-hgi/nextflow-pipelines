nextflow.preview.dsl=2

// list of cram from local disk, attached to a dummy study_id:
params.local_crams = "$baseDir/../../inputs/cram_list.csv" // set to "" if no local crams
// list of Irods study whose cram files must also be processed:
params.irods_studies = "$baseDir/../../inputs/study_ids_list.csv" // set to "" if no Irods crams

// required for mosdepth: genome ref fasta, capture region and study alias name of each study_id:
params.study_ref_capture = "$baseDir/../../inputs/study_ref_capture.csv"

params.dropqc = ""
include baton_study_id from '../modules/mosdepth/baton.nf' params(run: true, outdir: params.outdir)
include mosdepth_from_irods from '../modules/mosdepth/mosdepth_from_irods.nf' params(run:true, outdir: params.outdir,dropqc: params.dropqc)
include mosdepth_from_local from '../modules/mosdepth/mosdepth_from_local.nf' params(run:true, outdir: params.outdir,dropqc: params.dropqc)
include depth_plots_per_study from '../modules/mosdepth/depth_plots_per_study.nf' params(run:true, outdir: params.outdir)
include depth_plots from '../modules/mosdepth/depth_plots.nf' params(run:true, outdir: params.outdir)
include prep_capture_bed from '../modules/mosdepth/prep_capture_bed.nf' params(run:true, outdir: params.outdir)

workflow {
    // prepare bed file: remove comments and ^chr if required
    Channel.fromPath(params.study_ref_capture)
	.splitCsv(header: true, sep: ',')
	.map{row->tuple(row.study_id, file(row.capture), row.remove_capture_chr)}
	.set{to_prep_capture_bed}
    to_prep_capture_bed.view()
    prep_capture_bed(to_prep_capture_bed)
    
    Channel.fromPath(params.study_ref_capture)
	.splitCsv(header: true, sep: ',')
	.map{row->tuple(row.study_id, row.study_alias, file(row.ref))}
	.combine(prep_capture_bed.out.study_capturebed, by: 0)
	.set{ch_study_alias_ref_capture}
    ch_study_alias_ref_capture.view()
    
    //// process cram files from local_disk:
    ch_local_crams = params.local_crams == '' ? Channel.empty() :
	Channel.fromPath(params.local_crams)
	.splitCsv(header: true, sep: ',')
	.map{row->tuple(row.study_id, row.samplename, file(row.cram))}

    ch_local_crams
	.combine(ch_study_alias_ref_capture, by: 0)
	.set{to_mosdepth_local}
    // to_mosdepth_local.view()
    mosdepth_from_local(to_mosdepth_local.take(10))

    //// process cram files from irods:
    ch_irods_studies = params.irods_studies == '' ? Channel.empty() :
	Channel.fromPath(params.irods_studies)
	.splitCsv(header: true, sep: ',')
	.map{row-> row.study_id}
    
    // ch_irods_studies.view()
    baton_study_id(ch_irods_studies)
    baton_study_id.out.samples_noduplicates_tsv
	.map{a,b -> b}
	.splitCsv(header: true, sep: '\t')
	.map{row->tuple(row.sample, row.sample_supplier_name, row.study_id)}
	.map{a,b,c-> tuple(c,a)}
	.combine(ch_study_alias_ref_capture, by: 0)
	.set{to_mosdepth_irods}
    // to_mosdepth_irods.view()
    mosdepth_from_irods(to_mosdepth_irods.take(10))

    mosdepth_from_irods.out.study_samplename_regiondist
	.mix(mosdepth_from_local.out.study_samplename_regiondist)
	.set{ch_study_samplename_regiondist}
//	ch_study_samplename_regiondist
//	    .map{a,b,c-> tuple(a,c)}
//	    .groupTuple())

    depth_plots_per_study(
	ch_study_samplename_regiondist
	    .map{a,b,c-> tuple(a,c)}
	    .groupTuple())

    // depth_plots()
}
