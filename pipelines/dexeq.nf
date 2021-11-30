nextflow.preview.dsl=2

params.min_reads = 500   // used by crams_to_fastq_gz
//params.genome = 'GRCh38' // used by star aligner

params.fcextra = ""      // used by featurecounts
params.min_pct_aln  = 5 // used to filter STAR alignements, checking if align rate below threshold
params.singleend = false       // used by featurecounts
params.forward_stranded = false  // used by featurecounts

params.unstranded = false  // used by featurecounts
params.biotypes_header= "$baseDir/../assets/biotypes_header.txt" // used by featurecounts
params.mito_name = 'MT' // used by mapsummary
params.runtag = 'dddmouse' // HG_RNASeq of cellular DDD models pilot
params.ensembl_lib = "Ensembl 99 EnsDb" // used by tximport, must match used genome version
params.dropqc = ""
params.run_deseq2 = false
params.deseq2_tsv = "$baseDir/../../inputs/DESeq2.tsv"

params.run_star = true
def pick_aligner(String aligner) {
    return  aligner == 'star' || (!params.run_star && aligner == 'hisat2')
    ? true
    : false }

// params.pe_suffix_pattern  = '_{1,2}.fastq.gz'
// params.se_suffix    = '.fastq.gz'

//params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
params.star_index = "/lustre/scratch118/humgen/resources/conda/star_index"
Channel.fromPath(params.star_index)
    .ifEmpty { exit 1, "star index file not found: ${params.star_index}" }
    .set { ch_star_index}

// params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.gtf = "/lustre/scratch118/humgen/resources/ddd_mouse_genomes/with_cassette/Mus_musculus.GRCm38.99.gtf"
Channel.fromPath(params.gtf)
    .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
    .set { ch_gtf_star }

params.gff_dexseq = "/lustre/scratch118/humgen/resources/conda/DEXSeq.Mus_musculus.GRCm38.99.gff"
Channel.fromPath(params.gff_dexseq)
    .ifEmpty { exit 1, "GFF annotation file not found: ${params.gtf}" }
    .set { ch_gff_dexseq }

// params.salmon_index = params.genome ? params.genomes[ params.genome ].salmon ?: false : false
// params.salmon_index = "/lustre/scratch115/projects/interval_wgs/nextflow/salmon14_index/salmon"
params.salmon_index = "/lustre/scratch118/humgen/resources/conda/salmon_tindex"
Channel.fromPath(params.salmon_index)
    .ifEmpty { exit 1, "Salmon index dir not found: ${params.salmon_index}" }
    .set {ch_salmon_index}


Channel.fromPath(params.biotypes_header)
    .ifEmpty { exit 1, "biotypes header file not found: ${params.biotypes_header}" }
    .set { ch_biotypes_header }

// params.salmon_trans_gene = params.genome ? params.genomes[ params.genome ].salmon_trans_gene ?: false : false

//params.salmon_trans_gene = "/lustre/scratch115/projects/interval_wgs/nextflow/salmon14_index/trans_gene.txt"
//Channel.fromPath(params.salmon_trans_gene)
//    .ifEmpty { exit 1, "Salmon trans gene file not found: ${params.salmon_trans_gene}" }
//    .set {ch_salmon_trans_gene}


include iget_cram from '../modules/rna_seq/irods.nf' params(run:true, outdir: params.outdir,
								    dropqc: params.dropqc)
include crams_to_fastq_gz from '../modules/rna_seq/crams_to_fastq.nf' params(run:true, outdir: params.outdir,
								    min_reads: params.min_reads)
include fastqc from '../modules/rna_seq/fastqc.nf' params(run: true, outdir: params.outdir)

include salmon from '../modules/rna_seq/salmon.nf' params(run: true, outdir: params.outdir)
include merge_salmoncounts from '../modules/rna_seq/merge_salmoncounts.nf' params(run: true, outdir: params.outdir,
						   runtag : params.runtag)
include tximport from '../modules/rna_seq/tximport.nf' params(run: true, outdir: params.outdir,
						 ensembl_lib: params.ensembl_lib)

include star_2pass_basic from '../modules/rna_seq/star_2pass_basicmode.nf' params(run: true, outdir: params.outdir)
include star_transcriptomesam from '../modules/rna_seq/star_transcriptomesam.nf' params(run: true, outdir: params.outdir)
include dexseq_transcript_count from '../modules/rna_seq/dexseq_transcript_count.nf' params(run: true, outdir: params.outdir)

include star_2pass_1st_pass from '../modules/rna_seq/star_2pass_firstpass.nf' params(run: true, outdir: params.outdir)
include star_2pass_merge_junctions from '../modules/rna_seq/star_2pass_merge_junctions.nf' params(run: true, outdir: params.outdir)
include star_2pass_2nd_pass from '../modules/rna_seq/star_2pass_secondpass.nf' params(run: true, outdir: params.outdir)

include filter_star_aln_rate from '../modules/rna_seq/filter_star_aln_rate.nf' params(run: true,
									     min_pct_aln: params.min_pct_aln)
include leafcutter_bam2junc from '../modules/rna_seq/leafcutter_bam2junc.nf' params(run: true, outdir: params.outdir)
include leafcutter_clustering from '../modules/rna_seq/leafcutter_clustering.nf' params(run: true, outdir: params.outdir)
include featureCounts from '../modules/rna_seq/featurecounts.nf' params(run: true,outdir: params.outdir,
							       fcextra: params.fcextra,
							       singleend: params.singleend, 
							       forward_stranded: params.forward_stranded,
							       reverse_stranded: params.reverse_stranded,
							       unstranded: params.unstranded)
include samtools_index_idxstats from '../modules/rna_seq/samtools_index_idxstats.nf' params(run: true, outdir: params.outdir)
include mapsummary from '../modules/rna_seq/mapsummary.nf' params(run: true, outdir: params.outdir,
							 mito_name: params.mito_name)
include merge_featureCounts from '../modules/rna_seq/merge_featureCounts.nf' params(run: true, outdir: params.outdir,
									   runtag : params.runtag)
include multiqc from '../modules/rna_seq/multiqc.nf' params(run: true, outdir: params.outdir,
						   runtag : params.runtag)
include lostcause from '../modules/rna_seq/lostcause.nf' params(run: true, outdir: params.outdir,
						   runtag : params.runtag)
include baton_study_id from '../modules/rna_seq/baton.nf' params(run: true, outdir: params.outdir,
						   runtag : params.runtag)
include heatmap from '../modules/rna_seq/heatmap.nf' params(run: true, outdir: params.outdir,
						   runtag : params.runtag)
include deseq2 from '../modules/rna_seq/deseq2.nf' params(run: true, outdir: params.outdir, runtag: params.runtag)
include star_tabgenes_matrix from '../modules/rna_seq/star_tabgenes_matrix.nf' params(run: true, outdir: params.outdir, runtag: params.runtag)


workflow {

    to_dexseq = Channel.fromPath('/lustre/scratch114/projects/chromo_ndd/nextflow_rna_seq_dddmouse/results/star_pass2_basic/bams_starbasic.tsv')
	.splitCsv(header: true, sep: ',')
	.map{row->tuple(row.samplename, file(row.bam))}
	.map{a,b-> tuple(a,b)}
	//.toSortedList()
    
    dexseq_transcript_count(to_dexseq, ch_gff_dexseq.collect())

    // if(params.run_deseq2)
//	deseq2(salmon.out[0].collect(), Channel.fromPath(params.deseq2_tsv))
       
}
