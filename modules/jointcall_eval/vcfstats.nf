params.run = true

process vcfstats {
    //tag ""
    memory = '4G'
    time '240m'
    cpus 1
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    maxForks 12
    publishDir "${params.outdir}/${tag}/vcfstats/", mode: 'symlink', overwrite: true
    conda '/lustre/scratch118/humgen/resources/conda_envs/R3.6'

    when:
    params.run

    input: 
    set file(vcf), file(tbi)
    val(tag)

    output: 
    file("plots") //, emit: vcf_tbi

    script:
    """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/R3.6/bin:\$PATH
R_LIBS_SITE=/lustre/scratch118/humgen/resources/conda_envs/R3.6/lib/R/library
R_LIBS_USER=/lustre/scratch118/humgen/resources/conda_envs/R3.6/lib/R/library

mkdir plots

vcfstats --vcf $vcf \\
	--outdir plots \\
	--formula 'COUNT(1) ~ CONTIG' \\
	--title 'Number of variants on each chromosome'

vcfstats --vcf $vcf \\
	--outdir plots \\
	--formula 'COUNT(1, VARTYPE[snp]) ~ SUBST[A>T,A>G,A>C,T>A,T>G,T>C,G>A,G>T,G>C,C>A,C>T,C>G]' \\
	--title 'Number of substitutions of SNPs'

vcfstats --vcf $vcf \\
	--outdir plots \\
	--formula 'AAF ~ CONTIG' \\
	--title 'Allele frequency on each chromosome' 

vcfstats --vcf $vcf \\
	--outdir plots \\
	--formula 'COUNT(1, group=VARTYPE) ~ CHROM' \\
	--title 'Types of variants on each chromosome'

vcfstats --vcf $vcf \\
	--outdir plots \\
	--formula 'COUNT(1, group=VARTYPE) ~ CHROM[1]' \\
	--title 'Types of variants on each chromosome 1' \\
	--figtype pie

vcfstats --vcf $vcf \\
	--outdir plots \\
	--formula 'COUNT(1, group=GTTYPEs[HET,HOM_ALT]{0}) ~ CHROM' \\
        --title 'Mutant genotypes on each chromosome'
    """
}
// \\	--config examples/config.toml
