process ccs {
    publishDir "${params.outdir}", mode: 'copy'

    conda 'bioconda::pbccs=6.0.0'

    input:
    path bam

    output:
    path 'ccs.bam', emit: bam

    shell:
    """
    ccs ${bam} ccs.bam
    """
}