params.run = true 

process 'star_2pass_merge_junctions' {
    tag "STAR2_merge_junc"
    memory = '20 G'
    cpus 1
    time '100m'
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    
    publishDir "${params.outdir}/star2pass_merge_junc/", mode: 'copy', pattern: "SJ.filtered.tab"

  input:
    file(tabs) // from ch_tab_tabs_only.collect()

    when:
    params.run

    
  output:
    file "SJ.filtered.tab" //into ch_tab_filter_out

  script:

  """
  cat *.tab | awk '(\$5 > 0 && \$7 > 2 && \$6==0)' | cut -f1-6 | sort | uniq > SJ.filtered.tab
  """
}
