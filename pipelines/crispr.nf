nextflow.preview.dsl=2
params.runtag = 'walk66_subppol'
params.read2 = 'discard' // used by count_crispr_reads
params.min_reads = 500   // used by crams_to_fastq_gz

params.guide_libraries = "${baseDir}/../../guide_libraries/walkup66_check_subpool_library_cut.csv"
Channel.fromPath(params.guide_libraries)
    .set{ch_library_files}

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

//  //   1.A: from irods:
//    Channel.fromPath("${baseDir}/../../inputs/june_6007.csv")
//    	.splitCsv(header: true)
//    	.map { row -> tuple("${row.samplename}", "${row.batch}", "${row.sanger_sample_id}", "${row.study_id}") }
//    	.set{ch_to_iget}
//
//    iget_crams(ch_to_iget)
//    
//    iget_crams.out.iget_not_found
//	.map{ samplename, not_found_txt -> not_found_txt}
//	.collectFile(name: 'all_iget_not_found.txt', newLine: true, storeDir: "$params.outdir" )
//    
//    
//    crams_to_fastq_gz(iget_crams.out.spname_batch_cram) //.map{samplename,batch, crams,crais -> [samplename, batch, crams]})
//    crams_to_fastq_gz.out[0]
//	.map{samplename,batch, fastqs -> [samplename, batch, "1", fastqs]}
//	.set{ch_samplename_batch_starttrim_fastqs}


//
    // 1.B: or directly from fastq (if from basespace/lustre location rather than irods)
    Channel.fromPath("${baseDir}/../../inputs/walkup66_subpool.csv")
	.splitCsv(header: true)
	.map { row -> tuple("${row.samplename}", "${row.batch}", "${row.start_trim}", file("${row.fastq}")) }
	.set{ch_samplename_batch_fastqs}

    fastx_trimmer(ch_samplename_batch_starttrim_fastqs)

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
	.map{ samplename, fastqs -> tuple(samplename, fastqs, "June35.guide_library.csv", 0, 19, 0, 19) }
	.set{ch_samplename_fastq_library_match}
    
    count_crispr_reads(ch_samplename_fastq_library_match, ch_library_files.collect())

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

    
