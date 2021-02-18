#!/usr/bin/bash

# ##################################################
#Run this script: bash run_VCFCons.sh <sample>
#
#Expects to have the folllowing files:
#  <sample>.aligned.bam       <-- input to DeepVariant
#  <sample>.vcf               <-- output from DeepVariant
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

$SAMTOOLS index ${SAMPLE}.aligned.bam
$SAMTOOLS mpileup --min-BQ 1 -f $REF -s ${SAMPLE}.aligned.bam > ${SAMPLE}.aligned.bam.mpileup
$SAMTOOLS depth -q 0 -Q 0 ${SAMPLE}.aligned.bam > ${SAMPLE}.aligned.bam.depth

# IMPORTANT: for DeepVariant, we are turning off -q (QUAL) score filter in VCFCons.py by setting -q 0
# the scale of QUAL scores in DeepVariant is different
$PROG $REF $SAMPLE -c ${MIN_COVERAGE} -f ${MIN_ALT_FREQ} --vcf_type deepvariant -q 0 --input_depth ${SAMPLE}.aligned.bam.depth --input_vcf ${SAMPLE}.vcf

minimap2 -a $REF ${SAMPLE}.vcfcons.frag.fasta > ${SAMPLE}.vcfcons.frag.fasta.sam
$SAMTOOLS view -bS ${SAMPLE}.vcfcons.frag.fasta.sam > ${SAMPLE}.vcfcons.frag.fasta.bam
$SAMTOOLS sort ${SAMPLE}.vcfcons.frag.fasta.bam > ${SAMPLE}.vcfcons.frag.fasta.sorted.bam
$SAMTOOLS index ${SAMPLE}.vcfcons.frag.fasta.sorted.bam
