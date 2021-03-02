params.run = true

process 'leafcutter_clustering' {
    // cache 'deep'
    container "leafcutter"
    scratch '/tmp'
    stageInMode 'copy'
    stageOutMode 'rsync'
    queue 'long'
    memory ='100G'
    cpus 4
    time '2800m'
    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5
    
    publishDir "${params.outdir}/leafcutter/clustering", mode: 'symlink', pattern: "*.junc.clust.sorted.gz"
    publishDir "${params.outdir}/leafcutter/clustering", mode: 'symlink', pattern: "clust_*"
    publishDir "${params.outdir}/leafcutter/clustering", mode: 'symlink', pattern: "juncfiles.txt"
    publishDir "${params.outdir}/leafcutter/clustering", mode: 'symlink', pattern: "fofn_junctions_files.txt"

    input:
    file (junc_files) //from star_bam2junc.collect()

    when:
    params.run
    
  output:
    set file('clust_perind.counts.gz'), file('clust_perind_numers.counts.gz'), file('clust_pooled'),file('clust_refined'),file('clust_sortedlibs'), file('*.junc.clust.sorted.gz') //into star_clustering
    file("fofn_junctions_files.txt")
    // file "*.ReadsPerGene.out.tab" into ch_merge_starcounts

    script:
  """
  export PATH=/home/leafcutter/scripts:/home/leafcutter/clustering:\$PATH

  ls . | grep .junc\$ > fofn_junctions_files.txt
  leafcutter_cluster.py -j fofn_junctions_files.txt -m 50 -o clust -l 500000
  """
}

//  python /home/leafcutter/clustering/leafcutter_cluster.py -j ${juncfiles_txt} -m 50 -o clust -l 500000

process 'leafcutter_clustering_regtools' {
    // cache 'deep'
    container "leafcutter"
    scratch '/tmp'
    stageInMode 'copy'
    stageOutMode 'rsync'
    queue 'long'
    memory ='100G'
    cpus 4
    time '2800m'
    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5
    
    publishDir "${params.outdir}/leafcutter_regtools/clustering", mode: 'symlink', pattern: "*.junc.clust.sorted.gz"
    publishDir "${params.outdir}/leafcutter_regtools/clustering", mode: 'symlink', pattern: "clust_*"
    publishDir "${params.outdir}/leafcutter_regtools/clustering", mode: 'symlink', pattern: "juncfiles.txt"
    publishDir "${params.outdir}/leafcutter_regtools/clustering", mode: 'symlink', pattern: "fofn_junctions_files.txt"

    input:
    file (junc_files) //from star_bam2junc.collect()

    when:
    params.run
    
  output:
    set file('clust_perind.counts.gz'), file('clust_perind_numers.counts.gz'), file('clust_pooled'),file('clust_refined'),file('clust_sortedlibs'), file('*.junc.clust.sorted.gz') //into star_clustering
    file("fofn_junctions_files.txt")
    // file "*.ReadsPerGene.out.tab" into ch_merge_starcounts

    script:
  """
  export PATH=/home/leafcutter/scripts:/home/leafcutter/clustering:\$PATH

  ls . | grep .junc\$ > fofn_junctions_files.txt
  # pre regtools leafcutter leafcutter_cluster.py -j fofn_junctions_files.txt -m 50 -o clust -l 500000
  
  python2 /home/leafcutter/clustering/leafcutter_cluster_regtools.py -j fofn_junctions_files.txt --checkchrom -m 50 -r clust -o clust -l 500000  
  """
}
