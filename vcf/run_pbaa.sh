#!/usr/bin/bash

# ##################################################
#
# INPUT: a post-patient-demux BAM file (still has the amplicon primers)
# OUTPUT: pbaa VCF file <SAMPLE>.vcf
#
# ##################################################

INPUT_BAM=$1   # ex: <full_path>/demultiplexing_files/demultiplex.M13_bc1027_F--M13_bc1075_R.bam
SAMPLE=$2      # ex: LC00012
PRIMERS=/data/sarscov2_primers.fasta 
GENOME=/data/NC_045512.2.fasta
PBAA_REF=/data/sarscov2_guide_for_pbaa.fasta
PROG2TABLE=consensusVariants.py
PROG2VCF=pbaa2vcf.py


bamtools convert -format fastq -in sample.bam > sample.fastq
samtools faidx sample.fastq

# run pbaa
pbaa cluster --min-cluster-read-count 2 --trim-ends 0 $PBAA_REF --log-file pbaa.log sample.fastq $SAMPLE
# convert pbaa outcome to VCF
$PROG2TABLE -r $SAMPLE -p $SAMPLE --read_info ${SAMPLE}_read_info.txt --hifiSupport sample.fastq $GENOME ${SAMPLE}_passed_cluster_sequences.fasta
$PROG2VCF --passOnly -s barcode -o ${SAMPLE}.vcf ${SAMPLE}_alleles.csv ${SAMPLE}_variants.csv $GENOME
rm -rf ${SAMPLE}*.dot
