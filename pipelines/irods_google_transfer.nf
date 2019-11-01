#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                        A simple 1-step pipe to TRANSFER irods content to google 
========================================================================================
*/

def helpMessage() {
    log.info"""
    =========================================
     transfer from irods to google  pipeline v${version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf/transfer.nf --samples_and_irods_locs --google_bucket wsi_autozyg_part6_cram -profile farm4

    Mandatory arguments:
      -profile			Hardware config to use. farm3 / farm4 / docker
      --samples_and_irods_locs	File with sample names and irods locations
      --google_bucket		target google bucket 
      --outdir		 	directory for downloaded crams	
    """.stripIndent()
}


version = '1.0'

params.help = false
if (params.help){
    helpMessage()
    exit 0
}

params.samples_and_irods_locs = "/lustre/scratch118/humgen/hgi/projects/autozyg/Broad_gvcfs/Part6/irods_files_and_targets.csv"
params.google_bucket = "gs://wsi_autozyg_part6_cram"

if (!params.samples_and_irods_locs || !params.google_bucket) {
  exit 1, "Need both file containing samples_and_irods_locs and target google_bucket"
}

Channel
	.fromPath(params.samples_and_irods_locs)
	.splitCsv(header:true)
	.map{ row-> tuple(row.irods_file_location, row.target)}
	.set { ch_samples }

// Header log info
log.info "========================================="
log.info " transfer from irods to google  pipeline v${version} "
log.info "========================================="

def summary = [:]
summary['samples and irods locations']           = params.samples_and_irods_locs
summary['target google bucket']        = params.google_bucket

log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="



process 'copy_from_irods_and_copy_to_google' {
    executor 'local'
    maxForks 5
   
    publishDir "${params.outdir}/cram/", mode: 'move', pattern: "*.cram"
    publishDir "${params.outdir}/cram/", mode: 'move', pattern: "*.crai"

    input:
    set val(irods_file_location), val(target) from ch_samples

    println ">>> "+ch_samples+"<<<"
    
    output:
    set val(irods_file_location) into transferred 

    script:
    """
  
    # copy out of irods and copy to google
    # path looks like this: /seq/illumina/runs/30/30702/lane1/plex2/30702_1#2.cram 
    iget -v $irods_file_location $target
    /software/hgi/google-cloud-sdk/bin/gsutil cp $target ${params.google_bucket} 

    """
}
