nextflow.preview.dsl=2

params.runtag = 'deepvariant'
params.index_crams = false
params.run_deepvariant = true
params.outdir = "${baseDir}/../../outputs"
//params.ref_dir = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/"
//params.bed_dir = "/lustre/scratch118/humgen/resources/exome/Homo_sapiens/"

include sort_cram from '../modules/variant_calling/sort_cram.nf' params(run: true, outdir: params.outdir)
include markDuplicates from '../modules/variant_calling/markDuplicates.nf' params(run: true, outdir: params.outdir)
include coord_sort_cram from '../modules/variant_calling/coord_sort_cram.nf' params(run: true, outdir: params.outdir)
include deepvariant from '../modules/variant_calling/deepvariant.nf' params(run: true, outdir: params.outdir)

workflow {

    if (params.run_deepvariant) {
	

	ch_cram_file = Channel
		.fromPath("${baseDir}/../../inputs/cram_paths.txt")
		.splitText()
		.take(1)
		//.view()
     main:
	sort_cram(ch_cram_file)
	markDuplicates(sort_cram.out)
	coord_sort_cram(markDuplicates.out)
	deepvariant(coord_sort_cram.out)
     emit:
        my_data = deepvariant.out
        
    }
}
