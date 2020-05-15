#!/usr/bin/env python
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

def generate_cov_bed(count_filename, amplicon_bed):

    amplicon_info = {} #a amplicon name --> (0-based start, 1-based end)
    for line in open(amplicon_bed):
        ref, s0, e1, name = line.strip().split()
        amplicon_info[name] = int(s0), int(e1)

    cov = np.zeros(GENOME_SIZE, dtype=int)

    reader = DictReader(open(count_filename),delimiter='\t')
    for r in reader:
        n1, n2 = r['IdxFirstNamed'],r['IdxCombinedNamed']
        p1, p2 = n1.split('_')[0], n2.split('_')[0]
        assert p1 == p2  # lima must be called with --neighbor and with primer fastas that have agreeing prefix
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
    parser.add_argument("amplicon_bed", help="expected amplicon BED file")

    args = parser.parse_args()

    generate_cov_bed(args.lima_count, args.amplicon_bed)