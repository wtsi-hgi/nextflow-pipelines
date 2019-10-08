SELECT 
    manual_qc,
    study.id_study_lims,
    study.name,
    sample.common_name,
    sample.description,
    sample.name,
    sample.sanger_sample_id,
    sample.supplier_name,
    sample.donor_id,
    sample.gender,
    iseq_flowcell.last_updated,
    iseq_product_metrics.num_reads,
    iseq_product_metrics.id_run,
    iseq_product_metrics.position,
    iseq_product_metrics.tag_index
FROM
    iseq_flowcell
        JOIN
    sample ON iseq_flowcell.id_sample_tmp = sample.id_sample_tmp
        JOIN
    study ON iseq_flowcell.id_study_tmp = study.id_study_tmp
        JOIN
    iseq_product_metrics ON iseq_flowcell.id_iseq_flowcell_tmp = iseq_product_metrics.id_iseq_flowcell_tmp
WHERE
    study.id_study_lims IN (5591)
LIMIT 20000