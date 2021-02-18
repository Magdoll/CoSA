#!/usr/bin/bash

# ##################################################
#Run this script: bash run_VCFCons.sh <sample>
#
#Expects to have the folllowing files:
#  <sample>.bam
#  <sample>.vcf
#
#Will output:
#  <sample>.vcfcons.fasta
#  <sample>.vcfcons.frag.fasta
# ##################################################

SAMPLE=$1

PROG=VCFCons.py
REF=/pbi/dept/bifx/etseng/projects2020/SARS_Cov2/NC_045512.2.fasta
MIN_COVERAGE=4
MIN_ALT_FREQ=0.5
SAMTOOLS=/home/UNIXHOME/etseng/software_downloads/samtools-1.11/samtools

$SAMTOOLS index ${SAMPLE}.bam
$SAMTOOLS mpileup --min-BQ 1 -f $REF -s ${SAMPLE}.bam > ${SAMPLE}.bam.mpileup
$SAMTOOLS depth -q 0 -Q 0 ${SAMPLE}.bam > ${SAMPLE}.bam.depth

$PROG $REF $SAMPLE -c ${MIN_COVERAGE} -f ${MIN_ALT_FREQ} --vcf_type pbaa --input_depth ${SAMPLE}.depth --input_vcf ${SAMPLE}.vcf

minimap2 -a $REF ${SAMPLE}.vcfcons.frag.fasta > ${SAMPLE}.vcfcons.frag.fasta.sam
$SAMTOOLS view -bS ${SAMPLE}.vcfcons.frag.fasta.sam > ${SAMPLE}.vcfcons.frag.fasta.bam
$SAMTOOLS sort ${SAMPLE}.vcfcons.frag.fasta.bam > ${SAMPLE}.vcfcons.frag.fasta.sorted.bam
$SAMTOOLS index ${SAMPLE}.vcfcons.frag.fasta.sorted.bam
