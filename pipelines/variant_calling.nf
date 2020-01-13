nextflow.preview.dsl=2
params.runtag = 'test'
params.run_graphtyper_pipeline = true
params.run_graphtyper = true

Channel.fromPath("${baseDir}/../../inputs/bamlist.txt")
	.set{ch_bamlist_file}

Channel.fromPath("${baseDir}/../../inputs/graphtyper_pipeline_config.sh")
	.set{ch_graphtyper_pipeline_config}



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

workflow {

    if (params.run_graphtyper_pipeline) {
	graphtyper_pipeline(ch_bamlist_file, ch_graphtyper_pipeline_config)

	graphtyper_pipeline.out.commands_split
	    .splitText()
	    .take(4)
	    .set{ch_commands_split}

	ch_commands_split.view()

	if (params.run_graphtyper) {
	    graphtyper(ch_bamlist_file.collect(), ch_graphtyper_pipeline_config.collect(), ch_commands_split)
	}
    }
}
