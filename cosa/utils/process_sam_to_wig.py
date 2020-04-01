import os, sys
from csv import DictReader
import numpy as np
from cupcake.io.BioReaders import GMAPSAMReader

REF_LENGTH = 29903

def process_sam_to_wig(sam_filename, output_wig, cov_threshold=200, meta_info=None):
    cov = np.zeros(REF_LENGTH)
    reader = GMAPSAMReader(sam_filename, True)

    f_sam = open(sam_filename[:sam_filename.rfind('.')] + '.metainfo.sam', 'w')
    f_sam.write(reader.header)
    bad_count = 0
    for r in reader:
        tags = ''
        if len(r.segments) > 1:
            for e in r.segments: cov[e.start:e.end] += 1
            bad_count += 1
        tags += "\tsg:i:{0}".format(len(r.segments))  # sg: number of segments
        if meta_info is not None:
            seqid = r.qID.split('|')[0]
            if seqid in meta_info:
                tags += "\tst:A:{0}".format(meta_info[seqid]['Sequencing technology'][0]) # st: sequencing technology
            else:
                print("WARNING: Could not find {0} in metadata. Skipping.".format(seqid))
        f_sam.write(r.record_line + tags + '\n')
    f_sam.close()

    for i in range(len(cov)):
        if cov[i] < cov_threshold: cov[i] = 0

    f = open(output_wig, 'w')
    f.write("variableStep chrom=NC_045512v2 start=1")
    for i in range(len(cov)):
        f.write("{0} {1}\n".format(i+1, cov[i]))
    f.close()

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("sam_filename")
    parser.add_argument("output_wig")
    parser.add_argument("-c", "--coverage_threshold", type=int, default=200, help="Threshold coverage under which to show 0")
    parser.add_argument("-m", "--metadata_csv", help="Metadata CSV (optional)")
    args = parser.parse_args()

    if args.metadata_csv is not None:
        meta_info = dict((r['Accession ID'], r) for r in DictReader(open(args.metadata_csv), delimiter=','))

    process_sam_to_wig(args.sam_filename, args.output_wig, args.coverage_threshold, meta_info)