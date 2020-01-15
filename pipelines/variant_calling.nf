nextflow.preview.dsl=2
params.runtag = 'iwes'

params.use_interval_list = true
params.run_graphtyper_on_interval = false
params.index_crams = true
Channel.fromPath("${baseDir}/../../inputs/iwes_intervals.csv")
	.set{ch_iwes_intervals_csv}
Channel.fromPath("${baseDir}/../../inputs/graphtyper_pipeline_config_on_interval.sh")
	.set{ch_graphtyper_pipeline_config}

params.run_graphtyper_pipeline = false 
params.run_graphtyper = false

Channel.fromPath("${baseDir}/../../inputs/bqsr_crams_4070.txt")
	.set{ch_bamlist_file}

// https://gitlab.internal.sanger.ac.uk/hgi-projects/ibd-x10/cromwell/blob/master/ScatterMarkDup_BQSR_HC_inputs.json
// "ref_dict": "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/hs38DH.dict",
// "ref_fasta": "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/hs38DH.fa",
// "ref_fasta_index": "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/hs38DH.fa.fai",

//params.genome_fasta = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/hs38DH.fa"
//Channel.fromPath(params.ch_genome_fasta)
//    .set {ch_genome_fasta}
//params.genome_fasta_fai = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/hs38DH.fa.fai"
//Channel.fromPath(params.ch_genome_fasta_fai)
//    .set {ch_genome_fasta_fai}

include graphtyper_pipeline from '../modules/variant_calling/graphtyper_pipeline.nf' params(run: true, outdir: params.outdir)
include graphtyper from '../modules/variant_calling/graphtyper.nf' params(run: true, outdir: params.outdir)
include graphtyper_on_interval from '../modules/variant_calling/graphtyper_on_interval.nf' params(run: true, outdir: params.outdir)
include index_cram from '../modules/variant_calling/index_cram.nf' params(run: true, outdir: params.outdir)

workflow {

    if (params.run_graphtyper_pipeline) {
	graphtyper_pipeline(ch_bamlist_file, ch_graphtyper_pipeline_config)

	graphtyper_pipeline.out.commands_split
	    .splitText()
	    .map{a -> a.replaceAll(~/$\s/, "")}
	    .map{a -> tuple(a, a.replaceAll(~/:.*$/, "").replaceAll(~/^.*chr/, "chr"))}
	    .take(-1)
	    .set{ch_commands_split}

	if (params.run_graphtyper) {
	    graphtyper(ch_bamlist_file.collect(), ch_graphtyper_pipeline_config.collect(), ch_commands_split)
	}
    }

    if (params.use_interval_list) {

	if(params.index_crams) {
	    index_cram(ch_bamlist_file
		       .splitText().take(-1).map{a -> file(a.replaceAll(~/$\s/, ""))})
	}
	
	ch_iwes_intervals_csv
	    .splitCsv(header: true)
	    .map { row -> tuple(row.chr, row.start, row.end)}
	    .take(100)
	    .set{ch_chr_start_end}

	if (params.run_graphtyper_on_interval) {
	    graphtyper_on_interval(ch_bamlist_file.collect(), ch_graphtyper_pipeline_config.collect(), ch_chr_start_end)
	}
    }


    
}
