process bamtools_merge {
    publishDir "${params.outdir}", mode: 'copy'

    conda 'bioconda::lima=2.5.1'

    input:
    tuple val(sample), val(bam)

    output:
    path "${sample}.bam", emit: bam

    shell:
    """
    bamtools merge -out ${sample}.bam -in ${bam}
    """
}