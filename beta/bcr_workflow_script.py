#!/usr/bin/env python

"""

sample_metadata.txt example:

Sample,PrimerF,PrimerR
IgG,D702_F,D501_R
IgM,D702_F,D502_R
IgK,D702_F,D503_R
IgL,D702_F,D504_R
IgGsub,D705_F,D501_R


Command example:

python bcrtcr_process_post_lima.py output.D702_F--D501_R.bam.fastq IgG -n -w 50 > IgG.log
python bcrtcr_process_post_lima.py output.D702_F--D501_R.bam.fastq IgM -n -w 50 > IgM.log
"""
import os, sys
from csv import DictReader

VALID_REP_NAMES = {'IGG':50,'IGM':50,'IGK':50,'IGL':50,'IGA':50,'IGGSUB':120}
UMI_CUTOFFS = [2, 3, 5, 10]

def read_metadata(filename):
    reader = DictReader(open(filename), delimiter=',')
    if 'Sample' not in reader.fieldnames or \
        'PrimerF' not in reader.fieldnames or \
        'PrimerR' not in reader.fieldnames:
        print("metadata file {0} must have headers Sample,PrimerF,PrimerR. Abort!", file=sys.stderr)
        sys.exit(-1)
    d = {}  # (F,R) --> sample
    for r in reader:
        if r['Sample'].upper() not in VALID_REP_NAMES:
            print("{0} is not a valid sample name! Must be one of {1}".format(r['Sample'], list(VALID_REP_NAMES.keys())), file=sys.stderr)
            sys.exit(-1)
        d[r['PrimerF'],r['PrimerR']] = r['Sample']
    return d

def main(args):

    f_cmd = open(args.output_cmd, 'w')
    meta_info = read_metadata(args.metadata)
    good_files = []
    for (pF,pR),sample in meta_info.items():
        bam_file = "{0}.{1}--{2}.bam".format(args.lima_prefix, pF, pR)
        if not os.path.exists(bam_file):
            print("WARNING: expected {0} but not found. Ignoring.".format(bam_file), file=sys.stderr)
        else:
            f_cmd.write("bamtools convert -format fastq -in {0} > {0}.fastq\n".format(bam_file))
            good_files.append((bam_file, sample))

    for bam_file,sample in good_files:
        f_cmd.write("python bcrtcr_process_post_lima.py {0}.fastq {1} -n -w {2} > {1}.log\n".format(
            bam_file, sample, VALID_REP_NAMES[sample.upper()]))

    for bam_file,sample in good_files:
        determined_fq = sample + '_' + sample.upper()[:3] + '.determined.fq'
        for cutoff in UMI_CUTOFFS:
            f_cmd.write("python bcr_filter_by_umi_count.py {0} -c {1}\n".format(determined_fq, cutoff))

    f_cmd.close()


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-p", "--lima_prefix", help="Lima output prefix")
    parser.add_argument("-m", "--metadata", help="Sample metadata file")
    parser.add_argument("-o", "--output_cmd", help="Output command file")

    args = parser.parse_args()
    main(args)