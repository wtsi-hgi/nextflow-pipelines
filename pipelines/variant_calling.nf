nextflow.preview.dsl=2
params.runtag = 'test'
params.run_graphtyper = true


// https://gitlab.internal.sanger.ac.uk/hgi-projects/ibd-x10/cromwell/blob/master/ScatterMarkDup_BQSR_HC_inputs.json
// "ref_dict": "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/hs38DH.dict",
// "ref_fasta": "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/hs38DH.fa",
// "ref_fasta_index": "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/hs38DH.fa.fai",
params.genome_fasta = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/hs38DH.fa"
Channel.fromPath(params.ch_genome_fasta)
    .set {ch_genome_fasta}
params.genome_fasta_fai = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/hs38DH.fa.fai"
Channel.fromPath(params.ch_genome_fasta_fai)
    .set {ch_genome_fasta_fai}

Channel.fromPath("${baseDir}/../../inputs/bamlist.txt")
	.set{ch_bamlist_file}

Channel.fromPath("${baseDir}/../../inputs/bamlist.txt")
	.set{ch_bamlist_file}

include graphtyper from '../modules/variant_calling/graphtyper.nf' params(run: true, outdir: params.outdir)

workflow {

    if (params.run_graphtyper)
	graphtyper(ch_bamlist_file)
}
