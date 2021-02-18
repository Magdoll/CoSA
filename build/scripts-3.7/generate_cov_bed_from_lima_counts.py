#!/home/UNIXHOME/etseng/anacondaPy37/envs/anaCogentPy37/bin/python
import os, sys
import numpy as np
import math
import random
from csv import DictReader

GENOME_SIZE = 29904  # sars cov 2 genome size

"""
After running lima on amplicon data, generate a coverage BED file that can be used for IGV.
Input: 
   - output.lima.counts
   - expected amplicon bed file
"""

def generate_cov_bed(count_filename, amplicon_bed, is_bedpe=False):
    """
    :parram count_filename: count file that is [output].lima.counts
    :param amplicon_bed: BED or BEDPE file for the expected amplicons
    """

    amplicon_info = {} #a amplicon name --> (0-based start, 1-based end)
    for line in open(amplicon_bed):
        raw = line.strip().split()
        if is_bedpe:
            if len(raw)<7:
                print("Expected BEDPE format with 7 columns, saw {0} instead. Abort!".format(len(raw)))
                sys.exit(-1)
            ref, s0, e0, ref2, s1, e1, name = raw
            assert ref==ref2  # this should just be the sars cov 2 ref genome
        else:
            if len(raw)<4:
                print("Expected BED format with 4 columns, saw {0} instead. Abort!".format(len(raw)))
                sys.exit(-1)
            ref, s0, e1, name = raw
        amplicon_info[name] = int(s0), int(e1)

    cov = np.zeros(GENOME_SIZE, dtype=int)

    reader = DictReader(open(count_filename),delimiter='\t')
    for r in reader:
        n1, n2 = r['IdxFirstNamed'],r['IdxCombinedNamed']
        if is_bedpe:
            p1 = n1 + '--' + n2
        else:
            p1, p2 = n1.split('_')[0], n2.split('_')[0]
            assert p1 == p2  # lima must be called with --neighbor and with primer fastas that have agreeing prefix
        if p1 not in amplicon_info:
            print("WARNING: amplicon {0} not in BEDPE file {1}. Ignoring.".format(p1, amplicon_bed))
        else:
            s0, e1 = amplicon_info[p1]
            cov[s0:e1] +=  int(r['Counts'])

    # output coverage as a bedgraph
    # chrom start end value
    i = 0
    j = 1
    while i < GENOME_SIZE-1 and j < GENOME_SIZE:
        if cov[j]!=cov[i]:
            print("{0}\t{1}\t{2}\t{3}".format(ref, i, j, cov[i]), file=sys.stdout)
            i = j
            j = i + 1
        else:
            j += 1


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("lima_count", help="lima count file (ex: output.lima.counts")
    parser.add_argument("amplicon_bed", help="expected amplicon BED/BEDPE file. If BEDPE, use with --bedpe")
    parser.add_argument("--bedpe", default=False, action="store_true", help="BEDPE format instead of BED format")

    args = parser.parse_args()

    generate_cov_bed(args.lima_count, args.amplicon_bed, args.bedpe)