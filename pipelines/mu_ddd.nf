nextflow.preview.dsl = 2
// parent_read_trios.nf

params.samples1_file = "$baseDir/../../inputs/mosaic_IDs.csv"
params.trios_bam_vcf = "$baseDir/../../inputs/trios_bams_vcfs.txt"
params.r_script_parent_trio = "$baseDir/../../inputs/parent_reads_trio.R"

params.allDDD_irods_crams = "$baseDir/../../inputs/irods_cram_locations.txt"
// temporary, until cram files are actually moved from /lustre to Irods
// params.allDDD_lustre_crams = "$baseDir/../../inputs/allDDD_Source_CRAMs.csv"

params.run_parent_trio = true
include parent_trio from '../modules/mu_ddd/parent_trio.nf' params(run: true, outdir: params.outdir)

workflow {
    Channel.fromPath(params.trios_bam_vcf)
	.splitCsv(header: true, sep: '\t')
	.map{row->tuple(row.proband_sanger_id,
			row.father_sanger_id,
			row.mother_sanger_id)}
	.set{proband_father_mother_ids}

    Channel.fromPath(params.samples1_file)
	.splitCsv(header: true, sep: ',')
	.map{row->tuple(file(row.sample_file).getSimpleName().replaceAll(/_cands/, ""),
			file(row.sample_file))}
	.combine(proband_father_mother_ids, by: 0)
	.set{sample_file_father_mother}

    if (params.run_parent_trio) {
	parent_trio(sample_file_father_mother.take(2),	
		    Channel.fromPath(params.allDDD_irods_crams).collect(),
		    Channel.fromPath(params.r_script_parent_trio).collect())
    }
}

		    //Channel.fromPath(params.allDDD_lustre_crams).collect(),
