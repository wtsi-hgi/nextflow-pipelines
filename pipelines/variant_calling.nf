nextflow.preview.dsl=2
params.runtag = 'ibd_concat'

params.index_crams = false

// Channel.fromPath("${baseDir}/../../inputs/iwes_intervals_chr2.csv")
Channel.fromPath("${baseDir}/../../inputs/S04380110_Padded_merged.bed")
	.set{ch_intersect_bed}

// not on interval file-list
params.run_intersect_concat = true 
ch_vcfs_concat_prefix = "ibd"

// colnames: shard,vcf,tbi,coord,x1,x2,x3,batch
Channel.fromPath("${baseDir}/../../inputs/part4.csv")
	.set{ch_input_shards}

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

include sect_concat_vcfs from '../modules/variant_calling/sect_concat_vcfs.nf' params(run: true, outdir: params.outdir)

workflow {

    if (params.run_intersect_concat) {
	
	ch_input_shards
	    .splitCsv(header: true)
	    .take(1)
	    .map { row -> tuple(row.batch, file(row.vcf), row.coord)}
	    .map{a,b,c -> tuple(a,b.mklink("${baseDir}/../../results/vcfs/$batch/${c}.output.vcf.gz"))}
	    .set{ch_vcfs}
	
	ch_vcfs.view()
//	ch_input_shards
//	    .splitCsv(header: true)
//	    .map { row -> tuple(row.batch, file(row.tbi))}
//	    .set{ch_tbis}
//
//	ch_vcfs.mix(ch_tbis)
//	    .groupTuple()
//	    .take(1)
//	    .set{ch_by_50}
//	
//	sect_concat_vcfs(ch_by_50, ch_intersect_bed.collect())

//	run_intersect_concat(ch_batches, ch_intersect_bed)
    }
}
//
//    graphtyper_pipeline.out.commands_split
//	    .splitText()
//	    .map{a -> a.replaceAll(~/$\s/, "")}
//	    .map{a -> tuple(a, a.replaceAll(~/:.*$/, "").replaceAll(~/^.*chr/, "chr"))}
//	    .take(-1)
//	    .set{ch_commands_split}
//
//	if (params.run_graphtyper) {
//	    graphtyper(ch_bamlist_file.collect(), ch_graphtyper_pipeline_config.collect(), ch_commands_split)
//	}
//    }
//
//    if (params.use_interval_list) {
//
//	if(params.index_crams) {
//	    index_cram(ch_bamlist_file
//		       .splitText().take(-1).map{a -> file(a.replaceAll(~/$\s/, ""))})
//	}
//	
//	ch_iwes_intervals_csv
//	    .splitCsv(header: true)
//	    .map { row -> tuple(row.chr, row.start, row.end)}
//	    .take(-1)
//	    .filter { it[0] ==~ /chr[56]/} //.filter { it[1] ==~ /^[cC].*/}
//	    .set{ch_chr_start_end}
//
//	if (params.run_graphtyper_on_interval) {
//	    graphtyper_on_interval(ch_bamlist_file.collect(), ch_graphtyper_pipeline_config.collect(), ch_chr_start_end)
//	}
//    }
//
//    if (params.concat_vcfs) {
//	concat_vcfs(ch_vcfs_to_concat, ch_vcfs_concat_prefix)
//    }
//
//    
//}
//
