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


include crams_to_fastq_gz from '../modules/crams_to_fastq_anyflag.nf' params(run:true, outdir: params.outdir,
								    min_reads: params.min_reads)
include fastqc from '../modules/fastqc.nf' params(run: true, outdir: params.outdir)
include salmon from '../modules/salmon.nf' params(run: true, outdir: params.outdir)
include merge_salmoncounts from '../modules/merge_salmoncounts.nf' params(run: true, outdir: params.outdir,
						   runtag : params.runtag)
include tximport from '../modules/tximport.nf' params(run: true, outdir: params.outdir,
						 ensembl_lib: params.ensembl_lib)

include star_2pass_basic from '../modules/star_2pass_basicmode.nf' params(run: true, outdir: params.outdir)
include filter_star_aln_rate from '../modules/filter_star_aln_rate.nf' params(run: true,
									     min_pct_aln: params.min_pct_aln)
include leafcutter_bam2junc from '../modules/leafcutter_bam2junc.nf' params(run: true, outdir: params.outdir)
include leafcutter_clustering from '../modules/leafcutter_clustering.nf' params(run: true, outdir: params.outdir)
include featureCounts from '../modules/featurecounts.nf' params(run: true,outdir: params.outdir,
							       fcextra: params.fcextra,
							       singleend: params.singleend, 
							       forward_stranded: params.forward_stranded,
							       reverse_stranded: params.reverse_stranded,
							       unstranded: params.unstranded)
include samtools_index_idxstats from '../modules/samtools_index_idxstats.nf' params(run: true, outdir: params.outdir)
include mapsummary from '../modules/mapsummary.nf' params(run: true, outdir: params.outdir,
							 mito_name: params.mito_name)
include merge_featureCounts from '../modules/merge_featureCounts.nf' params(run: true, outdir: params.outdir,
									   runtag : params.runtag)
include multiqc from '../modules/multiqc.nf' params(run: true, outdir: params.outdir,
						   runtag : params.runtag)
include lostcause from '../modules/lostcause.nf' params(run: true, outdir: params.outdir,
						   runtag : params.runtag)
include iget from '../modules/irods.nf' params(run: true, outdir: params.outdir)

workflow {

    //Channel.fromPath('/lustre/scratch115/projects/bioaid/mercury_gn5/bioaid/inputs/to_iget4043.csv')
    Channel.fromPath('/lustre/scratch115/projects/bioaid/mercury_gn5/bioaid/inputs/to_iget.csv')
	.splitCsv(header: true)
	.map { row -> tuple("${row.samplename}", "${row.sample}", "${row.study_id}") }
	.set{ch_to_iget}

    iget(ch_to_iget)
    
    crams_to_fastq_gz(iget.out.map{samplename,crams,crais -> [samplename, crams]})
    
    crams_to_fastq_gz.out[0].map{a,b -> b}.set{fastq_to_publish}
    publish:
    fastq_to_publish to: '/lustre/scratch115/projects/interval_wgs/nextflow/walkups_x2_irods_17oct/',
	enabled: true, mode: 'copy', overwrite: true
    }
