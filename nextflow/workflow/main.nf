#!/usr/bin/env nextflow

// DSL2 is used
nextflow.enable.dsl=2

/*-----------------------------------------------------------------------------
  Pipeline Processes (includes)
-----------------------------------------------------------------------------*/

included { ccs                      } from '../modules/ccs/main'
included { lima                     } from '../modules/lima/main'
included { combine_demux_by_patient } from '../modules/combine_demux_by_patient/main'


/*-----------------------------------------------------------------------------
  Pipeline Parameters
-----------------------------------------------------------------------------*/

if (!params.subreads_bam) {
    exit 1, "[Pipeline Error] Missing parameter 'subreads_bam'\n"
}

if (!params.outdir) {
    exit 1, "[Pipeline Error] Missing parameter 'outdir'\n"
}

/*-----------------------------------------------------------------------------
  Main Workflow
-----------------------------------------------------------------------------*/

workflow {
    // [Channel] the input subreads BAM (eg. <movie>.subreads.bam or movie>..hifi_reads.bam)
    ch_in_subreads_bam = Channel.fromPath(params.subreads_bam)

    // [Process] call the consensus reads from the subreads
    ccs(ch_in_subreads_bam)

    // [Process] demultiplex the CCS reads
    lima(ccs.out.subreads_bam)

    // [Process] combine into per-patient data
    // TODO generate the tab-delimited data file
    // TODO run combine_demux_by_patient

    // [Process] trim amplicon primers
    // TODO: support --neigbhors or --different
    /*
    lima -j 30 --split-bam-named \
     --neighbors \
     --ccs \
     --min-score-lead 10 \
     --min-score 80
     ccs.bam amplicon_primers.fasta output.bam
     */

     // [Process] variant calling
     // TODO: support parameterizign variant callers
     // TODO: support bcftools
     // TODO: support DeepVariant
     // TODO: support pbaa

     // [Process] generate consensus sequence using VCFCons
     // TODO: samtools depth
     // TODO: run VCFCons

     // [Process] Assign lineages using Pangolin or Nextclade
     // TODO

     // TODO: other miscellaneous items
     // - downsample reads by amplicon
     // - generate per amplicon coverage BED file
     //

}