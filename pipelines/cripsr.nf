nextflow.preview.dsl=2
params.read2 = 'discard' // used by count_crispr_reads

// list irods study id and sample names:
//Channel.fromPath('/lustre/scratch115/projects/bioaid/mercury_gn5/bioaid/inputs/to_iget4043.csv')
params_sample_manifest_irods = '/lustre/scratch115/projects/bioaid/mercury_gn5/bioaid/inputs/to_iget.csv'
Channel.fromPath(params_sample_manifest_irods)
    .splitCsv(header: true)
    .map { row -> tuple("${row.samplename}", "${row.sample}", "${row.study_id}") }
    .set{ch_to_iget}

// collect library tables:
params.guide_libraries = '/lustre/scratch115/projects/bioaid/mercury_gn5/bioaid/inputs/*.guide_library.txt'
Channel.fromPath(params.guide_libraries)
    .set{ch_library_files}

// add guide library of each sample:
params.samplename_library = '/lustre/scratch115/projects/bioaid/mercury_gn5/bioaid/inputs/samplename_library.csv'
Channel.fromPath(params.samplename_library)
    .splitCsv(header: true)
    .map { row -> tuple("${row.samplename}", "${row.library}") }
    .set{ch_samplename_library}



include iget from '../modules/crispr/irods.nf' params(run: true, outdir: params.outdir)
include crams_to_fastq_gz from '../modules/crispr/crams_to_fastq_anyflag.nf' params(run:true, outdir: params.outdir,
								    min_reads: params.min_reads)
include count_crispr_reads from '../modules/crispr/count_crispr_reads.nf' params(run: true, outdir: params.outdir,
							 read2: params.read2)
include collate_crispr_counts from '../modules/crispr/collate_crispr_counts.nf' params(run: true, outdir: params.outdir)



workflow {

    iget(ch_to_iget)
    
    crams_to_fastq_gz(iget.out.map{samplename,crams,crais -> [samplename, crams]})
    
    crams_to_fastq_gz.out[0]
	.combine(ch_samplename_library, by: 0)
	.set{ch_samplename_fastq_library}
    
    count_crispr_reads(ch_samplename_fastq_library, ch_library_files.collect())

}

//    collate_crispr_counts(count_crispr_reads.out[0].collect())


    // publish output files
//    crams_to_fastq_gz.out[0].map{a,b -> b}.set{fastq_to_publish}
//    publish:
//    fastq_to_publish to: '/lustre/scratch115/projects/interval_wgs/nextflow/walkup_101/',
//	enabled: true, mode: 'copy', overwrite: true
//    collate_crispr_counts.out[0] to: '/lustre/scratch115/projects/interval_wgs/nextflow/walkup_101/',
//	enabled: true, mode: 'copy', overwrite: true
    
