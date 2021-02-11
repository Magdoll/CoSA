import os, sys, glob, vcf
import pdb
import phasing.io.SAMMPileUpReader as sp2

#PATTERN = sys.argv[1] # "*.vcfcons.vcf
#OUTFILE = sys.argv[2] # ../LCFall.bcftools_min_alt_freq_0.5.vcf_stats.txt

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("PATTERN")
parser.add_argument("OUTFILE")
parser.add_argument("--use_vcf_info", action="store_true", default=False)

args = parser.parse_args()

files = glob.glob(args.PATTERN)
f = open(args.OUTFILE, 'w')
f.write("sample\tpos\ttype\tlen\tdepth\talt_count\n")
for file in files:
    print(file)
    dirname, prefix = os.path.dirname(file), os.path.basename(file)
    prefix = prefix.split('.')[0]

    if not args.use_vcf_info:
        pileup_info = {}
        for r in sp2.MPileUpReader(os.path.join(dirname, prefix + '.bam.mpileup')): pileup_info[r.pos] = r
    for v in vcf.VCFReader(open(file)):
        _ref, _alt = str(v.REF), str(v.ALT[0]) # we'll ignore multi variants for now
        if v.is_snp:
            t = 'SUB'
        else:
             delta = len(_alt)-len(_ref)
             if delta==0: t = 'SUB'
             elif delta>0: t = 'INS'
             else: t = 'DEL'

        if args.use_vcf_info: # for pbaa-based VCF
            x = v.samples[0]
            read_cov = int(x.data.DP)
            GTs = x.data.GT.split('|')
            if len(GTs) == 1:
                if type(x.data.AD) is list:
                    print("ERROR: {0}:{1} does not have the matching number of genotypes and counts!".format(prefix, v.POS))
                    alt_count = int(x.data.DP)
                else:
                    alt_count = int(x.data.AD)
            else: # multiple genotypes
                if type(x.data.AD) is not list or len(x.data.AD)!=len(GTs):
                    print("ERROR: {0}:{1} does not have the matching number of genotypes and counts! Set alt_freq to 1 for now!".format(prefix, v.POS))
                    alt_count = int(x.data.DP)
                else:
                    alt_count = 0
                    # for now, we use the ALT0 count only, which is genotype 1
                    for gt,ad in zip(GTs, x.data.AD):
                        if gt=='1': alt_count += ad
        else:
            mrec = pileup_info[v.POS - 1]
            read_cov = mrec.cov
            if t=='SUB':
                alt_count = mrec.counts[_alt] + mrec.counts[_alt.lower()]
            elif t=='INS':
                #pdb.set_trace()
                tmp = '+' + str(delta)
                alt_count = mrec.counts[tmp+_alt[1:]] + mrec.counts[tmp+_alt[1:].lower()]
            else:
                #pdb.set_trace()
                tmp = str(delta)  # delta is already -<something>, don't need to add '-' sign
                alt_count = mrec.counts[tmp+_ref[1:]] + mrec.counts[tmp+_ref[1:].lower()]
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(prefix, v.POS, t, abs(delta) if t != 'SUB' else 1, read_cov, alt_count))
        
f.close()
