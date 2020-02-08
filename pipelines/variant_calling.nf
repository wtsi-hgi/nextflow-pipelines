nextflow.preview.dsl=2
params.runtag = 'iwes'
/// lfs quota -h -g interval_wes /lustre/scratch114

params.index_crams = false

// on interval file-list use:
params.use_interval_list = true
params.run_graphtyper_on_interval = true
// Channel.fromPath("${baseDir}/../../inputs/iwes_intervals_chr2.csv")
Channel.fromPath("${baseDir}/../../inputs/iwes_intervals_all.csv")
	.set{ch_iwes_intervals_csv}
Channel.fromPath("${baseDir}/../../inputs/graphtyper_scripts/*.sh")
	.set{ch_graphtyper_pipeline_config}


////////// how-to prep intervals from bed format:
//  cd /lustre/scratch114/projects/interval_wes/graphtyper_test/
//  echo "chr,start,end" > inputs/iwes_intervals_all.csv
//  cat /lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/exomes/S04380110_Padded_merged.bed | sed s'/\t/,/'g >> inputs/iwes_intervals_all.csv 
//////////

////////// how-to prep intervals from bed format:
// 14142 intervals for chr2:
//  cd /lustre/scratch114/projects/interval_wes/graphtyper_test/
//  echo "chr,start,end" > inputs/iwes_intervals_chr2.csv
//  cat /lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/intervals/Agilent_no_overlaps/chr2_S04380110_Padded+1_merged.bed | sed s'/\t/,/'g >> inputs/iwes_intervals_chr2.csv 
//////////

// not on interval file-list
params.run_graphtyper_pipeline = false 
params.run_graphtyper = false

// to concat vcfs after run done
params.concat_vcfs = false
// ch_vcfs_to_concat = "/lustre/scratch114/projects/interval_wes/graphtyper_test/results/graphtyper/results/chr1/"

// https://confluence.sanger.ac.uk/display/HGI/Interval+WES
ch_vcfs_to_concat = "/lustre/scratch118/humgen/hgi/projects/interval_wes/joint_calls/4070_QC_samples/output_vcf/"
ch_vcfs_concat_prefix = "gatk_chr"

//Channel.fromPath("${baseDir}/../../inputs/bqsr_crams_downsampled.txt")
Channel.fromPath("${baseDir}/../../inputs/bqsr_crams_4070.txt")
	.set{ch_bamlist_file}

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
include graphtyper_on_interval from '../modules/variant_calling/graphtyper_on_interval.nf' params(run: true, outdir: params.outdir)
include index_cram from '../modules/variant_calling/index_cram.nf' params(run: true, outdir: params.outdir)
include concat_vcfs from '../modules/variant_calling/concat_vcfs.nf' params(run: true, outdir: params.outdir)

workflow {

    if (params.run_graphtyper_pipeline) {
	graphtyper_pipeline(ch_bamlist_file, ch_graphtyper_pipeline_config)

	graphtyper_pipeline.out.commands_split
	    .splitText()
	    .map{a -> a.replaceAll(~/$\s/, "")}
	    .map{a -> tuple(a, a.replaceAll(~/:.*$/, "").replaceAll(~/^.*chr/, "chr"))}
	    .take(-1)
	    .set{ch_commands_split}

	if (params.run_graphtyper) {
	    graphtyper(ch_bamlist_file.collect(), ch_graphtyper_pipeline_config.collect(), ch_commands_split)
	}
    }

    if (params.use_interval_list) {

	if(params.index_crams) {
	    index_cram(ch_bamlist_file
		       .splitText().take(-1).map{a -> file(a.replaceAll(~/$\s/, ""))})
	}
	
	ch_iwes_intervals_csv
	    .splitCsv(header: true)
	    .map { row -> tuple(row.chr, row.start, row.end)}
	    .take(-1)
	    .filter { it[0] ==~ /chr[3456789]/} //.filter { it[1] ==~ /^[cC].*/}
	    .set{ch_chr_start_end}

	if (params.run_graphtyper_on_interval) {
	    graphtyper_on_interval(ch_bamlist_file.collect(), ch_graphtyper_pipeline_config.collect(), ch_chr_start_end)
	}
    }

    if (params.concat_vcfs) {
	concat_vcfs(ch_vcfs_to_concat, ch_vcfs_concat_prefix)
    }

    
}
