#!/usr/bin/bash

# ##################################################
#
# INPUT: a post-patient-demux BAM file (still has the amplicon primers)
# OUTPUT: pbaa VCF file <SAMPLE>.vcf
#
# ##################################################

SAMPLE=$1      # ex: LC00012

PRIMERS=/data/sarscov2_primers.fasta 
PBAA_REF=/data/sarscov2_guide_for_pbaa.fasta
GENOME=/data/NC_045512.2.fasta
PROG2TABLE=consensusVariants.py
PROG2VCF=pbaa2vcf.py

#bamtools convert -format fastq -in ${SAMPLE}.bam > ${SAMPLE}.fastq
samtools faidx ${SAMPLE}.fastq

# run pbaa
pbaa cluster --min-cluster-read-count 2 --trim-ends 0 $PBAA_REF --log-file ${SAMPLE}_pbaa.log ${SAMPLE}.fastq ${SAMPLE}_pbaa
# convert pbaa outcome to VCF
$PROG2TABLE -r ${SAMPLE}_pbaa -p ${SAMPLE}_pbaa --read_info ${SAMPLE}_pbaa_read_info.txt --hifiSupport ${SAMPLE}.fastq $GENOME ${SAMPLE}_pbaa_passed_cluster_sequences.fasta
$PROG2VCF --passOnly -s barcode -o ${SAMPLE}.vcf ${SAMPLE}_pbaa_alleles.csv ${SAMPLE}_pbaa_variants.csv $GENOME
rm -rf ${SAMPLE}*.dot
