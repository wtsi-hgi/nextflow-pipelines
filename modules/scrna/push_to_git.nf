params.run = true

process 'push_to_git' {
    tag "git-push $file_to_push"
    memory = '3G'
    time '120m'
    cpus 1
    maxForks 12
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 6
    cache false
    
    when:
    params.run 

    input:
    file(file_to_push)
    val(current_date)
    
    output:
    tuple val("${file_to_push.getName()}"), stdout, emit: stdout
    file("current_date.txt")
    env(WORK_DIR), emit: work_dir_to_remove

    script:
    """
WORK_DIR=\$PWD

sleep \$((1 + RANDOM % 100 + RANDOM % 100))
echo push minimal spreadsheet as tsv to gitlab

echo ${current_date} > current_date.txt

# remove empty lines
awk 'NF' ${file_to_push} > tmp
mv tmp ${file_to_push}

cp ${file_to_push} ../../../sync_status/
cd ../../.. && git pull && git add sync_status/${file_to_push}
git commit --allow-empty -m "auto_sync"
git push -u origin

cd \$WORK_DIR
echo done
   """
}

