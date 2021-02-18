import os, sys, glob, vcf
from cosa.vcf import VCFCons

#PATTERN = sys.argv[1] # "*.vcfcons.vcf
#OUTFILE = sys.argv[2] # ../LCFall.bcftools_min_alt_freq_0.5.vcf_stats.txt

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("PATTERN")
parser.add_argument("OUTFILE")
parser.add_argument("--vcf_type", choices=['pbaa', 'deepvariant', 'CLC', 'bcftools'], default=None, help="VCF format info, only used if --use_vcf_info is ON")

args = parser.parse_args()

files = glob.glob(args.PATTERN)
f = open(args.OUTFILE, 'w')
f.write("sample\tpos\ttype\tlen\tdepth\talt_count\n")
for file in files:
    print(file)
    dirname, prefix = os.path.dirname(file), os.path.basename(file)
    prefix = prefix.split('.')[0]

    for v in vcf.VCFReader(open(file)):
        _ref, _alt = str(v.REF), str(v.ALT[0]) # we'll ignore multi variants for now
        if v.is_snp:
            t = 'SUB'
        else:
             delta = len(_alt)-len(_ref)
             if delta==0: t = 'SUB'
             elif delta>0: t = 'INS'
             else: t = 'DEL'

        x = v.samples[0]

        if args.vcf_type == 'pbaa':
            read_cov = x.data.DP
            alt_count_dict = VCFCons.get_alt_count_pbaa(x, "{0}:{1}".format(prefix, v.POS))
            alt_count = alt_count_dict.most_common()[0][1]
        elif args.vcf_type == 'CLC':
            read_cov = x.data.DP
            alt_count_dict = VCFCons.get_alt_count_clc(len(v.ALT) + 1, x, "{0}:{1}".format(prefix, v.POS))
            alt_count = alt_count_dict.most_common()[0][1]
        elif args.vcf_type == 'bcftools':
            ##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
            read_cov = v.INFO['DP']
            alt_count = v.INFO['DP4'][2] + v.INFO['DP4'][3]
        else:
            read_cov = x.data.DP
            alt_count_dict = VCFCons.get_alt_count_std(len(v.ALT) + 1, x, "{0}:{1}".format(prefix, v.POS))
            alt_count = alt_count_dict.most_common()[0][1]

        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(prefix, v.POS, t, abs(delta) if t != 'SUB' else 1, read_cov, alt_count))
        
f.close()
