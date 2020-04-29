#!/usr/bin/env python
import os, sys
import numpy as np

GENOME_SIZE = 29904  # sars cov 2 genome size

def calc_coverage(cov_filename, amplicon_bed=None):
    cov = np.zeros(GENOME_SIZE, dtype=int)
    for line in open(cov_filename):
        ref, pos1, c = line.strip().split()
        cov[int(pos1)-1] = int(c)

    valid_ranges = []
    valid_cov = np.zeros(GENOME_SIZE, dtype=int)
    if amplicon_bed is not None:
        # bed format: chrom, start, end, name
        for line in open(amplicon_bed):
            ref, s0, e1, name = line.strip().split()
            s0, e1 = int(s0), int(e1)
            valid_ranges.append((s0,e1))
            valid_cov[s0:e1] += 1

    print("# avg coverage:", sum(cov)*1./GENOME_SIZE)
    n_miss_base  = sum(cov==0)
    print("missing bases:", n_miss_base, "(", round(n_miss_base*100./GENOME_SIZE,2), "%)")
    if amplicon_bed is not None:
        print("---accounting for expected amplicon coverage---")
        eff_len = sum(valid_cov>0)
        print("bases covered by amplicon design:", eff_len)
        print("# avg coverage:", sum(sum(cov[a:b]) for a,b in valid_ranges)/eff_len)
        n_miss_base =  sum(cov[i]==0 for i in valid_cov.nonzero()[0])
        print("missing bases:", n_miss_base, "(", round(n_miss_base*100./eff_len,2), "%)")

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("coverage_file", help="Coverage file in format [ref],[pos],[coverage]")
    parser.add_argument("-a", "--amplicon_bed", help="Expected amplicon BED file. If given, calculate effective coverage.")

    args = parser.parse_args()

    calc_coverage(args.coverage_file, args.amplicon_bed)
