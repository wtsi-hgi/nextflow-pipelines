nextflow.preview.dsl=2

params.min_reads = 500   // used by crams_to_fastq_gz
params.genome = 'GRCh38' // used by star aligner
params.fcextra = ""      // used by featurecounts
params.min_pct_aln  = 5 // used to filter STAR alignements, checking if align rate below threshold
params.singleend = false       // used by featurecounts
params.forward_stranded = false  // used by featurecounts
params.reverse_stranded = true  // used by featurecounts
params.unstranded = false  // used by featurecounts
params.biotypes_header= "$baseDir/assets/biotypes_header.txt" // used by featurecounts
params.mito_name = 'MT' // used by mapsummary
params.runtag = 'interval_basic' // used by mapsummary and multiqc
params.ensembl_lib = "Ensembl 91 EnsDb" // used by tximport, must match used genome version

params.run_star = true
def pick_aligner(String aligner) {
    return  aligner == 'star' || (!params.run_star && aligner == 'hisat2')
    ? true
    : false }

// params.pe_suffix_pattern  = '_{1,2}.fastq.gz'
// params.se_suffix    = '.fastq.gz'

params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
Channel.fromPath(params.star_index)
    .ifEmpty { exit 1, "star index file not found: ${params.star_index}" }
    .set { ch_star_index}
    
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
Channel.fromPath(params.gtf)
    .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
    .set { ch_gtf_star }

Channel.fromPath(params.biotypes_header)
    .ifEmpty { exit 1, "biotypes header file not found: ${params.biotypes_header}" }
    .set { ch_biotypes_header }

// params.salmon_index = params.genome ? params.genomes[ params.genome ].salmon ?: false : false
params.salmon_index = "/lustre/scratch115/projects/interval_wgs/nextflow/salmon14_index/salmon"
Channel.fromPath(params.salmon_index)
    .ifEmpty { exit 1, "Salmon index dir not found: ${params.salmon_index}" }
    .set {ch_salmon_index}

// params.salmon_trans_gene = params.genome ? params.genomes[ params.genome ].salmon_trans_gene ?: false : false
params.salmon_trans_gene = "/lustre/scratch115/projects/interval_wgs/nextflow/salmon14_index/trans_gene.txt"
Channel.fromPath(params.salmon_trans_gene)
    .ifEmpty { exit 1, "Salmon trans gene file not found: ${params.salmon_trans_gene}" }
    .set {ch_salmon_trans_gene}


include crams_to_fastq_gz from './modules/crams_to_fastq.nf' params(run:true, outdir: params.outdir,
								    min_reads: params.min_reads)
include fastqc from './modules/fastqc.nf' params(run: true, outdir: params.outdir)

include salmon from './modules/salmon.nf' params(run: true, outdir: params.outdir)
include merge_salmoncounts from './modules/merge_salmoncounts.nf' params(run: true, outdir: params.outdir,
						   runtag : params.runtag)
include tximport from './modules/tximport.nf' params(run: true, outdir: params.outdir,
						 ensembl_lib: params.ensembl_lib)

include star_2pass_basic from './modules/star_2pass_basicmode.nf' params(run: true, outdir: params.outdir)
include filter_star_aln_rate from './modules/filter_star_aln_rate.nf' params(run: true,
									     min_pct_aln: params.min_pct_aln)
include leafcutter_bam2junc from './modules/leafcutter_bam2junc.nf' params(run: true, outdir: params.outdir)
include leafcutter_clustering from './modules/leafcutter_clustering.nf' params(run: true, outdir: params.outdir)
include featureCounts from './modules/featurecounts.nf' params(run: true,outdir: params.outdir,
							       fcextra: params.fcextra,
							       singleend: params.singleend, 
							       forward_stranded: params.forward_stranded,
							       reverse_stranded: params.reverse_stranded,
							       unstranded: params.unstranded)
include samtools_index_idxstats from './modules/samtools_index_idxstats.nf' params(run: true, outdir: params.outdir)
include mapsummary from './modules/mapsummary.nf' params(run: true, outdir: params.outdir,
							 mito_name: params.mito_name)
include merge_featureCounts from './modules/merge_featureCounts.nf' params(run: true, outdir: params.outdir,
									   runtag : params.runtag)
include multiqc from './modules/multiqc.nf' params(run: true, outdir: params.outdir,
						   runtag : params.runtag)
include lostcause from './modules/lostcause.nf' params(run: true, outdir: params.outdir,
						   runtag : params.runtag)

workflow {

    Channel.fromPath('/lustre/scratch115/projects/interval_rna/inputs/*.cram').
    map{ it -> [ it.toString().replaceAll(~/.*\/(.*).cram/, "\$1"), it ] }.
    groupTuple(sort: true). //take(4).
    set{ch_cram_files}

    crams_to_fastq_gz(ch_cram_files)
    
    fastqc(crams_to_fastq_gz.out[0])

    salmon(crams_to_fastq_gz.out[0], ch_salmon_index.collect(), ch_salmon_trans_gene.collect())

    merge_salmoncounts(salmon.out[0].collect(), salmon.out[1].collect())

    tximport(salmon_out_0_1.map{it -> it.getName()}.collectFile(name: 'quant_sf_files.txt', sort: true, newLine: true), salmon_out_0_2)

    star_2pass_basic(crams_to_fastq_gz.out[0], ch_star_index.collect(), ch_gtf_star.collect())

    leafcutter_bam2junc(star_2pass_basic.out[0])

    leafcutter_clustering(leafcutter_bam2junc.out[0].collect())

    filter_star_aln_rate(star_2pass_basic.out[1].map{samplename,logfile,bamfile -> [samplename,logfile]}) // discard bam file, only STAR log required to filter
    
    filter_star_aln_rate.out.branch {
        filtered: it[1] == 'above_threshold'
        discarded: it[1] == 'below_threshold'}.set { star_2pass_basic_filter }
    
    star_2pass_basic_filter.filtered.combine(star_2pass_basic.out[1], by:0) //reattach bam file
	.map{samplename,filter,logfile,bamfile -> ["star", samplename, bamfile]} // discard log file and attach aligner name
	.set{star_2pass_basic_filtered} 
    
    samtools_index_idxstats(star_2pass_basic_filtered)
    
    mapsummary(samtools_index_idxstats.out)
    
    featureCounts(star_2pass_basic_filtered, ch_gtf_star.collect(), ch_biotypes_header.collect())

    merge_featureCounts(featureCounts.out[0].collect())
			//.collectFile(sort:true) { aligner, files -> [ aligner, files.collect{ it.toString() }.join('\n') + '\n' ] })
		       // .transpose()
		//	.groupTuple(sort: true)
		//	.map{ aligner, files -> [ aligner, files.collect{ it.toString() }.join('\n') + '\n' ] }
		//	.collectFile(sort:true) ) 

    crams_to_fastq_gz.out[1]
	.mix(star_2pass_basic_filter.discarded.map{samplename, filter -> [text: "${samplename}\tSTAR\tlowmapping\n"]})
	.set{ch_lostcause }
    lostcause(ch_lostcause.collectFile({ ['lostcause.txt', it.text]},sort:true))

    featureCounts.out[1]
	.filter{ pick_aligner(it[0]) }
	.map { it[1] }
	.set{ ch_multiqc_fc_aligner }

    featureCounts.out[2]
	.filter{ pick_aligner(it[0]) }
	.map{ it[1] }
	.set{ ch_multiqc_fcbiotype_aligner }
    
    multiqc(lostcause.out.collect().ifEmpty([]),
	    fastqc.out.collect().ifEmpty([]),
	    mapsummary.out.collect().ifEmpty([]),
	    ch_multiqc_fc_aligner.collect().ifEmpty([]),
	    ch_multiqc_fcbiotype_aligner.collect().ifEmpty([]),
	    star_2pass_basic.out[2].collect().ifEmpty([]),
	    salmon.out[2].collect().ifEmpty([]))
    
}
