process ccs {
    publishDir "${params.outdir}", mode: 'copy'

    conda 'bioconda::lima=2.0.0'

    input:
    path demux_bam
    path amplicon_primers_fasta

    output:
    path 'output.bam', emit: bams

    shell:
    def cores = 16
    if (task.cpus) {
        cores = (task.cpus as int)
    }
    def primers_arg = "--neighbors"
    if (params.lima_trim_amplicons_neighbors)      primers_arg = "--neighbors"
    else if (params.lima_trim_amplicons_different) primers_arg = "--different"
    else log.info "[lima trim amplicons] defaulting to use --neighbors"
    """
    lima \
        --num-threads ${cores} \
        -j 30 \
        --${primers_arg} \
        --ccs \
        --min-score-lead 10 \
        --min-score 80
        ${demux_bam} \
        ${amplicon_primers_fasta} \
        output.bam
    """
}