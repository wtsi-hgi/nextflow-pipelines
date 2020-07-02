params.run = true
params.runtag = 'runtag'

process collate_crispr_counts {
    tag "collate counts"
  //  conda '/software/hgi/installs/anaconda3/envs/nextflow20'
    // container "nfcore-rnaseq"
    // containerOptions = "--bind /lustre"
    publishDir "${params.outdir}/collate_counts", mode: 'symlink'
    memory = '20G'
    cpus 2
    time '300m'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    when:
    params.run

    input:
    set val(guide_library), file(samplename_counts_txt_files)

    output:
    set val(guide_library), file("*.count_matrix.txt")
    set val(guide_library), file("*.fofn_countsfiles.txt")

    shell:
    """
export PATH=/software/hgi/installs/anaconda3/envs/nextflow20/bin/:\$PATH

    echo count_file > fofn_files.txt
    ls . | grep .counts.txt\$ | grep -v genes.counts.txt >> fofn_files.txt

    echo samplename > fofn_samplenames.txt
    ls . | grep .counts.txt\$ | grep -v genes.counts.txt | sed s/.counts.txt// >> fofn_samplenames.txt

    paste -d ',' fofn_samplenames.txt fofn_files.txt > ${guide_library}.fofn_countsfiles.txt

    python3 ${workflow.projectDir}/../bin/crispr/collate_counts.py ${guide_library}.fofn_countsfiles.txt \"${guide_library}\"
    """
}
