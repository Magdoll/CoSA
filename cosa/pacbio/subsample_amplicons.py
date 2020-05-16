#!/usr/bin/env python
import os, sys
import math
import random
from csv import DictReader
# ex: primer name nCoV_2019_9_5p  nCoV_2019_9_alt2_3p

def downsample_lima_bam(count_filename, prefix, size=1000, valid_pairs_file=None):

    valid_combo = set()
    if valid_pairs_file is not None:
        for line in open(valid_pairs_file):
            n1, n2 = line.strip().split()
            valid_combo.add((n1,n2))  # don't care which order, add both
            valid_combo.add((n2,n1))

    good_counts = {}
    reader = DictReader(open(count_filename),delimiter='\t')
    for r in reader:
        n = "{0}.{1}--{2}".format(prefix,r['IdxFirstNamed'],r['IdxCombinedNamed'])
        #if r['IdxFirstNamed'].split('_')[2]==r['IdxCombinedNamed'].split('_')[2]:
        if valid_pairs_file is None or (r['IdxFirstNamed'],r['IdxCombinedNamed']) in valid_combo:
            good_counts[n] = int(r['Counts'])
        elif valid_pairs_file is not None:
            print("LOG: Ignoring invalid pairing {0}--{1} because not in whitelist".format(r['IdxFirstNamed'],r['IdxCombinedNamed']), file=sys.stderr)

    # sanity check all expected bams are there
    for k in good_counts:
        bam_file = k + '.bam'
        if not os.path.exists(bam_file):
            print("Expected lima output bam {0} but not found! Abort!".format(bam_file), file=sys.stderr)
            sys.exit(-1)

    for bam_prefix, counts in good_counts.items():
        frac = size / counts
        if frac >= 1:
            s = ''
        else:
            digit = -math.floor(math.log10(frac))
            s = "-s {seed}.{pad}{fraction}".format(seed=random.randint(1,10),
                                                   pad='0'*(digit-1),
                                                   fraction="{0:.0f}".format(int(frac*(10**digit))))
        #elif frac >= 0.01:  # frac 0.01-1.0
        #    s = "-s {seed}.0{fraction}".format(seed=random.randint(1,10),
        #                                      fraction=int(frac*100))
        #elif frac >= 0.001: # frac <= 0.01
        #    s = "-s {seed}.00{fraction}".format(seed=random.randint(1,10),
        #                                      fraction="{0:03.0f}".format(frac*100*1000))
        #print(size, counts, frac, s)
        #input()

        print("samtools view -b {s} {i}.bam | bamtools convert -format fastq > {i}.subsampled.fastq".format(\
                 s=s,
                 i=bam_prefix))

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("lima_prefix", help="lima output prefix")
    parser.add_argument("-s", "--subsample_size", type=int, default=1000, help="Number of reads to downsample to per amplicon (default: 1000)")
    parser.add_argument("--valid_pairs_file", help="(optional) file indicating list of valid primer pairs to subsample from")

    args = parser.parse_args()

    p = args.lima_prefix
    count_filename = p + '.lima.counts'

    if not os.path.exists(count_filename):
        print("Expected lima output file {0} but not found! Abort!".format(count_filename), file=sys.stderr)
        sys.exit(-1)

    if not os.path.exists(args.valid_pairs_file):
        print("File {0} not found! Abort!".format(args.valid_pairs+file), file=sys.stderr)
        sys.exit(-1)

    downsample_lima_bam(count_filename, args.lima_prefix, args.subsample_size, args.valid_pairs_file)
