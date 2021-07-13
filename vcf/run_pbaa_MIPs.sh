#!/usr/bin/bash

# ##################################################
#
# INPUT: a post-patient-demux BAM file (still has the amplicon primers)
# OUTPUT: pbaa VCF file <SAMPLE>.vcf
#
# ##################################################

SAMPLE=$1      # ex: LC00012

GENOME=NC_045512.2.fasta
CPUS=12
PROG2TABLE=~/GitHub/CoSA/vcf/consensusVariants.py

bamtools convert -format fastq -in ${SAMPLE}/alignment.bam > ${SAMPLE}/input.fastq
samtools faidx ${SAMPLE}/input.fastq

# run pbaa
pbaa cluster --log-level INFO --log-file ${SAMPLE}/input.pbaa.log --min-cluster-frequency 0.000001 --max-reads-per-guide 25000 --max-alignments-per-read 25000 -j $CPUS $GENOME ${SAMPLE}/input.fastq ${SAMPLE}/input.pbaa

# convert pbaa outcome to VCF
python $PROG2TABLE -r ${SAMPLE}/input.pbaa -p ${SAMPLE}/input.pbaa --hifiSupport ${SAMPLE}/input.fastq --read_info ${SAMPLE}/input.pbaa_read_info.txt --vcf --vcfMerge --vcfSampleCol runName --noCSV $GENOME ${SAMPLE}/input.pbaa_passed_cluster_sequences.fasta
