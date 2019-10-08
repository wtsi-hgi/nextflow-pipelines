Channel.fromPath('inputs/batch7_pheno.csv').
splitCsv(header: true).
map { row -> tuple("${row.samplename}", "${row.irods_path}") }.
map { a,b -> [a.toString(), b.toString()] }.
take( 2 ).
into{ch_csv_1;ch_csv_2;ch_csv_3;}

segment_number = 'batch7'
params.runtag = "iget_batch7"

params.run_irods = true

// map{ it -> [it, it.replaceAll(~/.*\/(.*).cram/, "\$1") ] }
process irods_iget {
    tag "iget ${segment_number} ${samplename}"
    memory '2G'
    cpus 1
    time '120m'
    errorStrategy 'terminate'
    //errorStrategy 'retry'
    //maxRetries 5

    maxForks 15 // maximum number of process instances that can be executed in parallel
    //publishDir "$params.outdir/input_crams/${segment_number}/", mode: 'move', pattern: "${samplename}.cram", overwrite: true
    publishDir "$params.outdir/input_crams/${segment_number}/", mode: 'symlink', pattern: "${samplename}.cram"
    publishDir "$params.outdir/input_crams/${segment_number}/", mode: 'symlink', pattern: "${samplename}.cram.crai"

    when:
    params.run_irods
    
    input: 
    set val(samplename), val(irod_file_path) from ch_csv_1

    output: 
    file("${samplename}.cram") into ch_cram
    file("${samplename}.cram.crai") into ch_cram_crai

    script:
    """ 
    iget -K ${irod_file_path}
    iget -K ${irod_file_path}.crai
    """
}
    //iget -K -n \$(ils -l ${irod_file_path} | awk '/green/ {if (\$4!=0) {print \$2; exit} }') ${irod_file_path} || iget -K ${irod_file_path}
    //iget -K -n \$(ils -l ${irod_file_path}.crai | awk '/green/ {if (\$4!=0) {print \$2; exit} }') ${irod_file_path}.crai || iget -K ${irod_file_path}.crai
