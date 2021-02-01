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

PROG=~/GitHub/CoSA/beta/VCFCons.py
REF=/pbi/dept/bifx/etseng/projects2020/SARS_Cov2/NC_045512.2.fasta
CDC_RENAME=/home/UNIXHOME/etseng/projects2021/SARS_Cov2/LabCorp/Jan2021_VCFcons_forCDC/LCFile_CDC_rename.txt
MIN_COVERAGE=4
MIN_ALT_FREQ=0.5

samtools index ${SAMPLE}.bam
samtools mpileup --min-BQ 1 -f $REF -s ${SAMPLE}.bam > ${SAMPLE}.bam.mpileup

python $PROG $REF $SAMPLE -c ${MIN_COVERAGE} -f ${MIN_ALT_FREQ} -s ${CDC_RENAME}

minimap2 -a $REF ${SAMPLE}.vcfcons.frag.fasta > ${SAMPLE}.vcfcons.frag.fasta.sam
samtools view -bS ${SAMPLE}.vcfcons.frag.fasta.sam > ${SAMPLE}.vcfcons.frag.fasta.bam
samtools sort ${SAMPLE}.vcfcons.frag.fasta.bam > ${SAMPLE}.vcfcons.frag.fasta.sorted.bam
samtools index ${SAMPLE}.vcfcons.frag.fasta.sorted.bam
