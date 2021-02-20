#!/usr/bin/bash

# ##################################################
#Run this script: bash run_VCFCons.sh <sample>
#
#Expects to have the folllowing files:
#  <sample>.aligned.bam            <-- input reads aligned to NC_045512.2.fasta
#  <sample>.vcf                    <-- output from pbaa + conversion to VCF
#
#Will output:
#  <sample>.vcfcons.fasta
#  <sample>.vcfcons.frag.fasta
# ##################################################

SAMPLE=$1

PROG=VCFCons.py
REF=/data/SARS_Cov2/NC_045512.2.fasta
MIN_COVERAGE=4
MIN_ALT_FREQ=0.5
SAMTOOLS=/software/samtools-1.11/samtools

$SAMTOOLS index ${SAMPLE}.aligned.bam
$SAMTOOLS depth -q 0 -Q 0 ${SAMPLE}.aligned.bam > ${SAMPLE}.aligned.bam.depth

$PROG $REF $SAMPLE -c ${MIN_COVERAGE} -f ${MIN_ALT_FREQ} --vcf_type pbaa --input_depth ${SAMPLE}.aligned.bam.depth --input_vcf ${SAMPLE}.vcf

minimap2 -a $REF ${SAMPLE}.vcfcons.frag.fasta > ${SAMPLE}.vcfcons.frag.fasta.sam
$SAMTOOLS view -bS ${SAMPLE}.vcfcons.frag.fasta.sam > ${SAMPLE}.vcfcons.frag.fasta.bam
$SAMTOOLS sort ${SAMPLE}.vcfcons.frag.fasta.bam > ${SAMPLE}.vcfcons.frag.fasta.sorted.bam
$SAMTOOLS index ${SAMPLE}.vcfcons.frag.fasta.sorted.bam
