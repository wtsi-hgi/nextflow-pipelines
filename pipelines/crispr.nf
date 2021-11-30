nextflow.preview.dsl=2
params.runtag = 'walkup103_v1.1.0'
params.read2 = 'discard' // used by count_crispr_reads
params.min_reads = 500   // used by crams_to_fastq_gz

// collect library tables:
// params.guide_libraries = "${baseDir}/../../guide_libraries/*.guide_library.csv"
params.guide_libraries = "${baseDir}/../../guide_libraries/tim_7nov.csv"
//params.guide_libraries = "${baseDir}/../../guide_libraries/tim_contamination.csv"
Channel.fromPath(params.guide_libraries)
    .set{ch_library_files}

// add guide library of each sample:
// params.samplename_library = "${baseDir}/../../inputs/walkup101_libraries.csv"
params.samplename_library = "${baseDir}/../../inputs/walkup103_libraries.csv"
Channel.fromPath(params.samplename_library)
    .splitCsv(header: true)
    .map { row -> tuple("${row.samplename}", "${row.library}", "${row.includeG}") }
    .set{ch_samplename_library}



include iget_crams from '../modules/crispr/irods.nf' params(run: true, outdir: params.outdir)
include crams_to_fastq_gz from '../modules/crispr/crams_to_fastq_anyflag.nf' params(run:true, outdir: params.outdir,
								    min_reads: params.min_reads)
include fastx_trimmer from '../modules/crispr/trim_fastq.nf' params(run: true, outdir: params.outdir)
include merge_fastq_batches from '../modules/crispr/merge_fastq_batches.nf' params(run:true, outdir: params.outdir)
include count_crispr_reads from '../modules/crispr/count_crispr_reads.nf' params(run: true, outdir: params.outdir,
							 read2: params.read2)
include collate_crispr_counts from '../modules/crispr/collate_crispr_counts.nf' params(run: true, outdir: params.outdir)
include fastqc from '../modules/crispr/fastqc.nf' params(run: true, outdir: params.outdir)
include multiqc from '../modules/crispr/multiqc.nf' params(run: true, outdir: params.outdir,
						   runtag : params.runtag)

workflow {

    // 1.A: from irods:
    //Channel.fromPath('/lustre/scratch115/projects/bioaid/mercury_gn5/bioaid/inputs/to_iget.csv')
    //	.splitCsv(header: true)
    //	.map { row -> tuple("${row.samplename}", "${row.batch}", "${row.sample}", "${row.study_id}") }
    //	.set{ch_to_iget}
    // iget_crams(ch_to_iget)
    // crams_to_fastq_gz(iget.out.map{samplename,batch, crams,crais -> [samplename, batch, crams]})
    // crams_to_fastq_gz.out[0].set{ch_samplename_batch_fastqs}

    // 1.B: or directly from fastq (if from basespace/lustre location rather than irods)
    // Channel.fromPath("${baseDir}/../../inputs/walkup101_fastqs.csv")
    Channel.fromPath("${baseDir}/../../inputs/walkup103_fastqs.csv")
	.splitCsv(header: true)
	.map { row -> tuple("${row.samplename}", "${row.batch}", "${row.start_trim}", file("${row.fastq}")) }
	.set{ch_samplename_batch_fastqs}

    fastx_trimmer(ch_samplename_batch_fastqs)

    fastx_trimmer.out 
	.groupTuple(by: 0, sort: true)
	.map{ samplename, batchs, fastqs -> tuple( groupKey(samplename, batchs.size()), batchs, fastqs ) }
	.set{ch_samplename_fastqs_to_merge}

    merge_fastq_batches(ch_samplename_fastqs_to_merge)

    fastqc(fastx_trimmer.out 
	    .map{ samplename, batch, fastq -> tuple( samplename, fastq ) }
	    .mix(merge_fastq_batches.out[0]))

    multiqc(fastqc.out.collect())

    merge_fastq_batches.out[0]
	.combine(ch_samplename_library, by: 0)
	.set{ch_samplename_fastq_library_includeG}
    
    count_crispr_reads(ch_samplename_fastq_library_includeG, ch_library_files.collect())

    collate_crispr_counts(
	count_crispr_reads.out[0]
	    .map{ lib_csv,counts -> [ lib_csv.replaceAll(~/.csv/, ""), counts ] }
	    .transpose()
	    .groupTuple(by: 0, sort: true)
	    .mix(count_crispr_reads.out[0].
		 map{lib,counts -> counts}.collect().map{a -> tuple("all_libs", a)})
    )
    
    count_crispr_reads.out[1].collectFile(name: 'mapping_percent.txt', newLine: true,
    					  storeDir: "${params.outdir}/", sort: true)

    // publish output files
//    crams_to_fastq_gz.out[0].map{a,b -> b}.set{fastq_to_publish}
//    publish:
//    fastq_to_publish to: '/lustre/scratch115/projects/interval_wgs/nextflow/walkup_101/',
//	enabled: true, mode: 'copy', overwrite: true
//    collate_crispr_counts.out[0] to: '/lustre/scratch115/projects/interval_wgs/nextflow/walkup_101/',
//	enabled: true, mode: 'copy', overwrite: true

}

    
