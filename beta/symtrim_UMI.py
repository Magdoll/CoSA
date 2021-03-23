#!/usr/bin/env python3
import os, sys
from Bio.Seq import Seq
from csv import DictReader, DictWriter
import pysam
import pdb



def clip_out(bam_filename, output_prefix, umi_len=5, extra_len=25, min_insert_len=300):
    """
    :param bam_filename: BAM of post-LIMA (primer-trimmed) CCS sequences

    M13-UMI-arm-insert-arm-UMI-M13
    (M13 should already be removed)
    """


    FIELDS = ['id', 'UMI1', 'UMI2', 'arm1', 'arm2', 'insert_len']

    f1 = open(output_prefix + '.trimmed.csv', 'w')
    writer1 = DictWriter(f1, FIELDS, delimiter=',')
    writer1.writeheader()

    reader = pysam.AlignmentFile(bam_filename, 'rb', check_sq=False)
    f2 = pysam.AlignmentFile(output_prefix+'.trimmed.bam', 'wb', header=reader.header)

    umi_extra_len = umi_len + extra_len

    for r in reader:
        d = r.to_dict()

        if len(d['seq'])-umi_extra_len*2 < min_insert_len:
            continue

        umi1 = d['seq'][:umi_len]
        umi2 = d['seq'][-umi_len:]
        arm1 = d['seq'][umi_len:umi_extra_len]
        arm2 = d['seq'][-umi_extra_len:-umi_len]

        rec = {'id': r.qname,
               'UMI1': umi1,
               'UMI2': umi2,
               'arm1': arm1,
               'arm2': arm2,
               'insert_len': len(d['seq'])-(umi_extra_len*2)}
        writer1.writerow(rec)

        d['seq'] = d['seq'][umi_extra_len:-umi_extra_len]
        d['qual'] = d['qual'][umi_extra_len:-umi_extra_len]
        new_tags = []
        for tag in d['tags']:
            if tag.startswith('zs:B'): # defunct CCS tag, don't use
                pass
            elif tag.startswith('dq:i:') or tag.startswith('iq:i:') or tag.startswith('sq:i:'):
                tag = tag[umi_extra_len:-umi_extra_len]
                new_tags.append(tag)
            else:
                new_tags.append(tag)

        new_tags.append('rg:Z:' + umi2)
        d['tags'] = new_tags
        x = pysam.AlignedSegment.from_dict(d, r.header)
        f2.write(x)

    f1.close()
    f2.close()
    print("Output written to: {0}.trimmed.csv|bam".format(output_prefix))


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("bam_filename", help="CCS BAM with cDNA primer removed (post LIMA)")
    parser.add_argument("output_prefix", help="Output prefix")
    parser.add_argument("-u", "--umi_len", type=int, default=5, help="Length of UMI (default: 5)")
    parser.add_argument("-e", "--extra_len", type=int, default=25, help="Length of arm sequence (default: 25)")
    parser.add_argument("--min_insert_len", type=int, default=300, help="Minimum insert length (default: 300)")


    args = parser.parse_args()

    if args.extra_len < 0:
        print("extra_len can't be a negative number!", file=sys.stderr)
        sys.exit(-1)
    if args.umi_len < 0:
        print("umi_len can't be a negative number!", file=sys.stderr)
        sys.exit(-1)

    clip_out(args.bam_filename,
             args.output_prefix,
             args.umi_len,
             args.extra_len,
             args.min_insert_len)