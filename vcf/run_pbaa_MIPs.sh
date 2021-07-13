#!/usr/bin/bash

# ##################################################
#
# INPUT: a post-patient-demux BAM file (still has the amplicon primers)
# OUTPUT: pbaa VCF file <SAMPLE>.vcf
#
# ##################################################

SAMPLE=$1      # ex: LC00012

GENOME=/home/UNIXHOME/etseng/projects2021/LabCorp_2021_COVID/July_more_samples_pangolin_share/10278_PBT5066_100004/NC_045512.2.fasta
CPUS=12
PROG2TABLE=consensusVariants.py

samtools faidx ${SAMPLE}/input.fastq

# run pbaa
pbaa cluster --log-level INFO --log-file ${SAMPLE}/input.pbaa.log --min-cluster-frequency 0.000001 --max-reads-per-guide 25000 --max-alignments-per-read 25000 -j $CPUS $GENOME ${SAMPLE}/input.fastq ${SAMPLE}/input.pbaa

# convert pbaa outcome to VCF
consensusVariants.py -r ${SAMPLE}/input.pbaa -p ${SAMPLE}/input.pbaa --hifiSupport ${SAMPLE}/input.fastq --read_info ${SAMPLE}/input.pbaa_read_info.txt --vcf --vcfMerge --vcfSampleCol runName --noCSV $GENOME ${SAMPLE}/input.pbaa_passed_cluster_sequences.fasta
