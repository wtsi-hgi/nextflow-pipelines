nextflow.preview.dsl=2
params.runtag = 'test'
params.run_graphtyper = true


params.genome_fasta = "$baseDir/../.fasta"
Channel.fromPath(params.ch_genome_fasta)
    .ifEmpty { exit 1, "genome fasta missing: ${params.genome_fasta}" }
    .set {ch_genome_fasta}


include graphtyper from '../modules/variant_calling/graphtyper.nf' params(run: true, outdir: params.outdir)

workflow {

    Channel.fromPath("${baseDir}/../../inputs/bamlist.txt")
	.set{ch_bamlist_file}

    if (params.run_graphtyper)
	graphtyper(ch_bamlist_file)
}
