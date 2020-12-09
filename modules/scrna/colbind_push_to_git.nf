params.run = true

process 'colbind_push_to_git' {
    tag "git-push samples_sync_statuts.tsv"
    memory = '3G'
    time '120m'
    cpus 1
    maxForks 12
    errorStrategy 'ignore'
    maxRetries 0
    publishDir "${params.outdir}/", mode: 'copy', pattern: "samples_sync_status.tsv"

    when:
    params.run 

    input:
    file(minimal_spreadsheet_tsv)
    file(sync_status_tsv)
    
    output:
    tuple val("samples_sync_status.tsv"), stdout, emit: stdout

    script:
    """
echo push minimal spreadsheet sync status as tsv to gitlab

echo remove empty lines
awk 'NF' ${sync_status_tsv} > ${sync_status_tsv}_rm_empty_lines
awk 'NF' ${minimal_spreadsheet_tsv} > ${minimal_spreadsheet_tsv}_rm_empty_lines

echo col bind
paste   | sed 's/\\t/\\0\\t/g' | column -s \$'\\t' -t > tmp

echo col bind
cat ${sync_status_tsv}_rm_empty_lines | sort > file1
cat ${minimal_spreadsheet_tsv}_rm_empty_lines | sort > file2
join -1 1 -2 1 -e0 -a1 -o'0,2.2,1.2,1.3' file2 file1 > tmp

cat tmp | grep sanger_sample_id > samples_sync_status.tsv
cat tmp | grep -v sanger_sample_id >> samples_sync_status.tsv

echo push to git 
cp samples_sync_status.tsv ../../../sync_status/
cd ../../.. && git pull && git add sync_status/samples_sync_status.tsv
git commit --allow-empty -m "auto_sync"
git push -u origin
echo done
   """
}


