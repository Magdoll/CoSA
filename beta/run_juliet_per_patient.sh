

PRIMER='/home/UNIXHOME/etseng/projects2020/SARS_Cov2/LabCorp/eden.primers.fasta'
GENOME='/home/UNIXHOME/etseng/projects2020/SARS_Cov2/LabCorp/NC_045512.2.fasta'
JULIET_CONFIG='/home/UNIXHOME/etseng/projects2020/SARS_Cov2/LabCorp/sarscov2.json'
CPUS=30

DIR=$1

pwd=$PWD


#cd $DIR
#lima -j $CPUS --split-bam-named --neighbors --ccs --min-score-lead -1 --peek-guess --window-size-multi 1.1 ccs.bam $PRIMER output.bam
#subsample_amplicons.py output --subsample_size 1000 > cmd
#bash cmd
#cat output*subsampled.fastq > subsampled.ccs.Q20.fastq
#rm -rf output*.xml
#rm -rf output*.pbi
pbmm2 align --sort --preset HiFi $GENOME subsampled.ccs.Q20.fastq subsampled.mapped.bam
juliet -c $JULIET_CONFIG subsampled.mapped.bam subsampled.minperc10.juliet.html subsampled.minperc10.juliet.json --min-perc 10
juliet_json_to_vcf.py subsampled.minperc10.juliet.json $JULIET_CONFIG subsampled.minperc10.juliet.vcf

# ------------------------------------------------------------
# consensus generation #1 
# using bcftools
# caveat: dropout regions will be using ref bases, no "N"s
# ------------------------------------------------------------
bgzip subsampled.minperc10.juliet.vcf
bcftools index subsampled.minperc10.juliet.vcf.gz
cat $GENOME|bcftools consensus subsampled.minperc10.juliet.vcf.gz > subsampled.minperc10.juliet.vcf_consensus.fasta

# ------------------------------------------------------------
# consensus generation #2
# using racon
# ------------------------------------------------------------
samtools view -h subsampled.mapped.bam > subsampled.mapped.sam
racon subsampled.ccs.Q20.fastq subsampled.mapped.sam $GENOME --no-trimming --no-trimming -q 20 -e 0.01 -t $CPUS > subsampled.minperc10.juliet.racon_consensus.fasta


cd $pwd
