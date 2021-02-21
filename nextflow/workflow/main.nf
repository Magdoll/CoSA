#!/usr/bin/env nextflow

// DSL2 is used
nextflow.enable.dsl=2

/*-----------------------------------------------------------------------------
  Pipeline Processes (includes)
-----------------------------------------------------------------------------*/

include { ccs                      } from '../modules/ccs/main'
include { lima                     } from '../modules/lima/main'
include { combine_demux_by_patient } from '../modules/combine_demux_by_patient/main'
include { bamtools_merge           } from '../modules/bamtools/merge/main'

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

    // [Channel] the sample metadata.  Should contain:
    //           "Sample", "BarcodeF", "BarcodeFName", "BarcodeR", "BarcodeRName"
    ch_in_metadata = Channel
        .fromPath(params.sample_metadata)
        .splitCsv(header: true)

    // [Channel] the barcode FASTA file, one entry per input barcode
    ch_in_barcodes_fasta = ch_in_metadata
        .map { row ->
            ">${row.BarcodeFName}\n${row.BarcodeF}\n>${row.BarcodeRName}\n${row.BarcodeR}"
        }
        .collectFile(name: 'barcodes.fasta', newLine: true)

    // [Process] call the consensus reads from the subreads
    ccs(ch_in_subreads_bam)

    // [Process] demultiplex the CCS reads
    lima(ccs.out.subreads_bam, ch_in_barcodes_fasta)

    // [Channel] build tuples of barcode key and BAM file
    ch_lima_bams = demux.out.bams.flatten().map { bam ->
        ["${bam.baseName}".substring("demux.".length()), bam]
    }

    // [[Channel] transform the metadata to tuples of barcode key and sample
    ch_in_lima_outputs = ch_in_metadata
        .map { row ->
            ["${row.BarcodeFName}--${row.BarcodeRName}", row.Sample]
        }

    // [[Channel] gather all BAMs for the same sample
    ch_bams_by_sample = ch_in_lima_outputs
        .join(ch_lima_bams)                        // join by barcode name key
        .map { key, sample, bam -> [sample, bam] } // discard the key
        .groupTuple()

    // [Process] combine into per-patient data
    bamtools_merge(ch_bams_by_sample)

    // [Process] trim amplicon primers
    // TODO: support --neigbhors or --different

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