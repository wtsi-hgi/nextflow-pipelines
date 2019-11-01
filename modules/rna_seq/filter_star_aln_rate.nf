params.run = true
params.min_pct_aln  = 5 // used to filter STAR alignements output if rate below threshold

process filter_star_aln_rate {
    tag "filter_star_aln_rate ${samplename}"

    //container 'nfcore-rnaseq' 
    errorStrategy 'retry'
    maxRetries 3
    time '30m'
    cpus 1
    memory '2G'

    //publishDir "${params.outdir}/fastq12/", mode: 'copy'
    
    when:
    params.run
    
    input:
    set val(samplename), file(log_final_out)
    output: 
    set val(samplename), stdout
    script:

    """
#!/usr/bin/env python3
import re

with open(\"${log_final_out}\",\"r\") as f:
    for line in f:
        if re.search(\"Uniquely mapped reads %\",line):
                if float(re.findall(\"\\d+\\.\\d+\", line)[0]) > float(${params.min_pct_aln}):
                    print('above_threshold', end='')
                else:
                    print('below_threshold', end='')
    """
}
