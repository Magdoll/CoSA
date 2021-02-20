#!/home/UNIXHOME/etseng/anacondaPy37/envs/anaCogentPy37/bin/python

import os, sys, subprocess
from csv  import DictReader
from collections import defaultdict

"""
Given lima output (after demux of, say, M13 barcodes) 
And sample metadata in tab-delimited format:

    BarcodeName Sample
    M13_bc1001_F--M13_bc1057_R  sample1
    M13_bc1002_F--M13_bc1057_R  sample1
    M13_bc1003_F--M13_bc1057_R  sample2
    M13_bc1004_F--M13_bc1057_R  sample2
    
This script will create, for each sample a directory and combine the lima output into a single XML or bam file.
"""

def combine_demux(lima_prefix, metadata_filename, fofn_only=False, force=False):

    reader = DictReader(open(metadata_filename), delimiter='\t')
    if 'BarcodeName' not in reader.fieldnames or 'Sample' not in reader.fieldnames:
        print("ERROR: expected fields BarcodeName and Sample in metadata file {0}!".format(metadata_filename), file=sys.stderr)
        sys.exit(-1)
    sample_info = dict((r['BarcodeName'], r['Sample']) for r in reader)
    files_per_sample = defaultdict(lambda: [])

    for bname,sname in sample_info.items():
        outfile = "{p}.{b}.bam".format(p=lima_prefix, b=bname)
        if not os.path.exists(outfile):
            print("WARNING: expected file {0} but not found. Ignoring...".format(outfile), file=sys.stderr)
        else:
            files_per_sample[sname].append(os.path.abspath(outfile))

    for sname in files_per_sample:
        if os.path.exists(sname) and not force:
            print("ERROR: directory {0} already exists.  Use --force to overwrite.", file=sys.stderr)
            sys.exit(-1)

    for sname, files in files_per_sample.items():
        print("Creating directory {0}...".format(sname),  file=sys.stdout)
        if not os.path.exists(sname):
            os.makedirs(sname)
        with open(os.path.join(sname, 'ccs.fofn'), 'w') as f:
            f.write("\n".join(files) + '\n')
            print("Writing fofn", f.name, file=sys.stdout)
        if not fofn_only:
            cmd = "bamtools merge -in " + " -in ".join(files) + " -out " + os.path.join(sname, 'ccs.bam')
            if subprocess.check_call(cmd, shell=True)==-1:
                print("ERROR running cmd:", cmd)
                sys.exit(-1)
            print("Writing merged bam", os.path.join(sname, 'ccs.bam'), file=sys.stdout)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("lima_prefix",  help="Lima output prefix after demux (ex: output)")
    parser.add_argument("sample_info", help="Barcode-Sample metadata")
    parser.add_argument("-f", "--force", default=False, action="store_true", help="Force overwriting existing directories. (default: off)")
    parser.add_argument("--fofn_only", default=False, action="store_true", help="Write out FOFN only(default: off)")

    args = parser.parse_args()
    combine_demux(args.lima_prefix, args.sample_info, args.fofn_only, args.force)





