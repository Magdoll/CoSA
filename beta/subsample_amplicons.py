import os, sys
from csv import DictReader
# ex: primer name nCoV_2019_9_5p  nCoV_2019_9_alt2_3p

def downsample_lima_bam(count_filename, size=1000):
    good_counts = {}
    reader = DictReader(open(count_filename),delimiter='\t')
    for r in reader:
        n = "output.{0}--{1}".format(r['IdxFirstNamed'],r['IdxCombinedNamed'])
        #if r['IdxFirstNamed'].split('_')[2]==r['IdxCombinedNamed'].split('_')[2]:
        good_counts[n] = int(r['Counts'])

    # sanity check all expected bams are there
    for k in good_counts:
        bam_file = k + '.bam'
        if not os.path.exists(bam_file):
            print("Expected lima output bam {0} but not found! Abort!".format(bam_file), file=sys.stderr)
            sys.exit(-1)

    for bam_prefix, counts in good_counts.items():
        a = int(min(1., size/counts)*100)
        s = ''
        if a < 100: s = "-s {seed}.{fraction}".format(seed=0, fraction=a)
        print("samtools view -b {s} {i}.bam | bamtools convert -format fastq > {i}.subsampled.fastq".format(\
                 s=s,
                 i=bam_prefix))

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("lima_prefix", help="lima output prefix")
    parser.add_argument("-s", "--subsample_size", type=int, default=1000, help="Number of reads to downsample to per amplicon (default: 1000)")

    args = parser.parse_args()

    p = args.lima_prefix
    count_filename = p + '.lima.counts'

    if not os.path.exists(count_filename):
        print("Expected lima output file {0} but not found! Abort!".format(count_filename), file=sys.stderr)
        sys.exit(-1)

    downsample_lima_bam(count_filename, args.subsample_size)