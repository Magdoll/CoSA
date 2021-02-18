#!/usr/bin/bash

# ##################################################
#Run this script: bash run_VCFCons.sh <sample>
#
#Expects to have the folllowing files:
#  <sample>.aligned.bam
#
# you can generate <sample>.aligned.bam using
#  pbmm2 align --sort --preset HiFi per_patient.primer_trimed.bam > sample.aligned.bam
#
#Will output:
#  <sample>.vcfcons.vcf
#  <sample>.vcfcons.fasta
#  <sample>.vcfcons.frag.fasta
# ##################################################

SAMPLE=$1

PROG=VCFCons.py
REF=/data/SARS_Cov2/NC_045512.2.fasta
MIN_COVERAGE=4
MIN_ALT_FREQ=0.5
SAMTOOLS=/software/samtools-1.11/samtools


samtools index ${SAMPLE}.aligned.bam
samtools mpileup --min-BQ 1 -f $REF -s ${SAMPLE}.aligned.bam > ${SAMPLE}.aligned.bam.mpileup
$SAMTOOLS depth -q 0 -Q 0 ${SAMPLE}.aligned.bam > ${SAMPLE}.aligned.bam.depth

# call variants using bcftools
bcftools mpileup -f $REF ${SAMPLE}.aligned.bam | bcftools call -mv -Ov -o ${SAMPLE}.vcf
# filter variants & generate consensus using VCFCons
$PROG $REF $SAMPLE -c ${MIN_COVERAGE} -f ${MIN_ALT_FREQ} --input_depth ${SAMPLE}.bam.depth --input_vcf ${SAMPLE}.vcf --vcf_type bcftools

minimap2 -a $REF ${SAMPLE}.vcfcons.frag.fasta > ${SAMPLE}.vcfcons.frag.fasta.sam
samtools view -bS ${SAMPLE}.vcfcons.frag.fasta.sam > ${SAMPLE}.vcfcons.frag.fasta.bam
samtools sort ${SAMPLE}.vcfcons.frag.fasta.bam > ${SAMPLE}.vcfcons.frag.fasta.sorted.bam
samtools index ${SAMPLE}.vcfcons.frag.fasta.sorted.bam
