process ccs {
    publishDir "${params.outdir}", mode: 'copy'

    conda 'bioconda::lima=2.0.0'

    input:
    path ccs_bam

    output:
    path 'demux.*', emit: bam

    shell:
    def cores = 16
    if (task.cpus) {
        cores = (task.cpus as int)
    }
    """
    lima --num-threads ${cores} \
        --split-bam-named \
        --different \
        --ccs \
        --min-score-lead 10 \
        --min-score 80 \
        ${ccs_bam} \
        demux.bam
    """
}