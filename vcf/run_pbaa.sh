#!/usr/bin/bash

SAMPLE=$1

REF=/home/UNIXHOME/etseng/projects2020/SARS_Cov2/LabCorp/Jan2021_VCFcons_forCDC/NC_045512.guide_for_pbaa.fasta
REF2=/home/UNIXHOME/etseng/projects2020/SARS_Cov2/LabCorp/Jan2021_VCFcons_forCDC/NC_045512.fa
PROG2TABLE=~/GitHub/CoSA/vcf/consensusVariants.py
PROG2VCF=~/GitHub/CoSA/vcf/pbaa2vcf.py

bamtools convert -format fastq -in ${SAMPLE}.bam > ${SAMPLE}.fastq
samtools faidx ${SAMPLE}.fastq
pbaa cluster --min-cluster-read-count 2 --trim-ends 0 --log-level DEBUG $REF --log-file ${SAMPLE}.log ${SAMPLE}.fastq $SAMPLE
python $PROG2TABLE -r $SAMPLE -p $SAMPLE --read_info ${SAMPLE}_read_info.txt --hifiSupport ${SAMPLE}.fastq $REF2 ${SAMPLE}_passed_cluster_sequences.fasta
python $PROG2VCF --passOnly -s barcode -o ${SAMPLE}.vcf ${SAMPLE}_alleles.csv ${SAMPLE}_variants.csv $REF2
rm -rf ${SAMPLE}*.dot
