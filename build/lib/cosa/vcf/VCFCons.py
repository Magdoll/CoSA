#!/usr/bin/env python3
__version__ = '8.1.0'
#import pdb
import os, sys
from collections import Counter
from csv import DictWriter
from Bio import SeqIO
import vcf

VARIANT_FIELDS = ['Pos', 'Type', 'Length', 'Depth', 'AltCount']

def make_seq_from_list(seqlist, start0, end1):
    seq = ''
    for i in range(start0, end1):
        if i in seqlist:
            seq += seqlist[i]
    return seq

def get_alt_count_std(num_gt, x, name):
    if num_gt != len(x.data.AD):
        print("ERROR: {0} does not have the matching number of genotypes and counts!".format(name))
        return Counter({'1': int(x.data.DP)})
    else:
        alt_count_dict = Counter()
        for i in range(1,num_gt):
            alt_count_dict[str(i)] += int(x.data.AD[i])
        return alt_count_dict

def get_alt_count_clc(num_gt, x, name):
    if num_gt != len(x.data.CLCAD2):
        print("ERROR: {0} does not have the matching number of genotypes and counts!".format(name))
        return Counter({'1':int(x.data.DP)})
    else:
        alt_count_dict = Counter()
        for i in range(1,num_gt):
            alt_count_dict[str(i)] += int(x.data.CLCAD2[i])
        return alt_count_dict

def get_alt_count_pbaa(num_gt, x, name):
    if num_gt == 1:
        if type(x.data.AD) is list:
            print("ERROR: {0} does not have the matching number of genotypes and counts!".format(name))
            return Counter({'1': int(x.data.DP)})
        else:
            return Counter({'1':int(x.data.AD)})
    else: # multiple genotypes
        if type(x.data.AD) is not list or len(x.data.AD)!=num_gt:
            print("ERROR: {0} does not have the matching number of genotypes and counts!".format(name))
            return Counter({'1': int(x.data.DP)})
        else:
            alt_count_dict = Counter()
            for i in range(1,num_gt):
                alt_count_dict[str(i)] += int(x.data.AD[i])
            return alt_count_dict

def genVCFcons(ref_fasta, depth_file, vcf_input, prefix, newid,
               min_coverage=4, min_alt_freq=0.5, min_qual=100,
               vcf_type=None):
    """
    :param ref_fasta: should be the Wuhan reference
    :param depth_file: <sample>.bam.depth of per base coverage
    :param vcf_input: VCF of where the variants are
    :param prefix: output prefix
    :param min_coverage: below this coverage bases will be 'N'
    :param min_alt_freq: below this ALT frequency bases will use the reference instead
    :param vcf_type: choices are pbaa, CLC, deepvariant (standard)
    :return:
    """
    output_fasta = prefix + '.vcfcons.fasta'
    output_frag_fasta = prefix + '.vcfcons.frag.fasta'
    output_info = prefix + '.vcfcons.info.csv'

    ref = next(SeqIO.parse(open(ref_fasta),'fasta'))
    refseq = str(ref.seq)
    refseq = list(refseq)

    # NOTE: we are using samtools depth to get the per base coverage
    depth_per_pos = {} # 0-based POS --> read depth
    for line in open(depth_file):
        chrom, pos1, count = line.strip().split()
        depth_per_pos[int(pos1)-1] = int(count)

    newseqlist = dict(zip(range(len(refseq)), list(refseq)))
    for pos0 in range(len(refseq)):
        if pos0 not in depth_per_pos or depth_per_pos[pos0] < min_coverage: newseqlist[pos0] = 'N'
    # make sure begin/ends are "N"s
    for pos0 in range(min(depth_per_pos)): newseqlist[pos0] = 'N'
    for pos0 in range(max(depth_per_pos),len(refseq)): newseqlist[pos0] = 'N'

    # now add in the variants
    tally_types = Counter() # SUB/INS/DEL --> count
    vcf_reader = vcf.Reader(open(vcf_input))
    vcf_writer = vcf.Writer(open(prefix+'.vcfcons.vcf', 'w'), vcf_reader)
    f_variant = open(prefix+'.vcfcons.variants.csv', 'w')
    variant_writer = DictWriter(f_variant, fieldnames=VARIANT_FIELDS, delimiter='\t')
    variant_writer.writeheader()
    for v in vcf_reader:
        # deepvariant has this weird record of RefCalls, ignore them
        if vcf_type == 'deepvariant' and v.FILTER == ['RefCall']: continue
        x = v.samples[0]
        # DeepVariant is unphased, can be 0/1, 1/1, etc...
        # pbaa is ?????
        if vcf_type == 'pbaa':
            total_cov = x.data.DP
            alt_count_dict = get_alt_count_pbaa(len(v.ALT)+1, x, "{0}:{1}".format(prefix, v.POS))
            alt_index, alt_count = alt_count_dict.most_common()[0]
        elif vcf_type == 'CLC':
            total_cov = x.data.DP
            alt_count_dict = get_alt_count_clc(len(v.ALT)+1, x, "{0}:{1}".format(prefix, v.POS))
            alt_index, alt_count = alt_count_dict.most_common()[0]
        elif vcf_type == 'bcftools':
            ##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
            total_cov = v.INFO['DP']
            alt_count = v.INFO['DP4'][2] + v.INFO['DP4'][3]
            alt_index = 1
        else:
            total_cov = x.data.DP
            alt_count_dict = get_alt_count_std(len(v.ALT)+1, x, "{0}:{1}".format(prefix, v.POS))
            alt_index, alt_count = alt_count_dict.most_common()[0]


        # alt_index is '1' for ALT0, '2' for ALT1...etc, so we have to do int(alt_index)-1 to get the genotype from v.ALT
        _ref, _alt = str(v.REF), str(v.ALT[int(alt_index)-1])
        if len(v.ALT)>1:
            print("WARNING: more than 1 alt type for {0}! Using just the ALT {1} cuz most abundant.".format(prefix, _alt))

        alt_freq = alt_count * 1. / total_cov
        _altlen = len(_alt)
        _reflen = len(_ref)
        delta = _altlen - _reflen
        if delta==0: t = 'SUB'
        elif delta>0: t = 'INS'
        else: t = 'DEL'

        if total_cov < min_coverage:
            print("INFO: For {0}: Ignore variant {1}:{2}->{3} because total cov is {4}.".format(prefix, v.POS, _ref, _alt, total_cov))
        elif alt_freq < min_alt_freq:
            print("INFO: For {0}: Ignore variant {1}:{2}->{3} because alt freq is {4}.".format(prefix, v.POS, _ref, _alt, alt_freq))
        elif v.QUAL is not None and v.QUAL < min_qual:
            print("INFO: For {0}: Ignore variant {1}:{2}->{3} because qual is {4}.".format(prefix, v.POS, _ref, _alt, v.QUAL))
        else:
            if v.QUAL is None:
                print("WARNING: QUAL field is empty for {0}:{1}. Ignoring QUAL filter.".format(prefix, v.POS))
            vcf_writer.write_record(v)
            variant_writer.writerow({'Pos': v.POS,
                                     'Type': t,
                                     'Length': abs(delta) if t != 'SUB' else 1,
                                     'Depth': total_cov,
                                     'AltCount': alt_count})
            tally_types[t] += 1
            if t=='SUB':
                # remember there could be consecutive subs
                for cur in range(_altlen):
                    newseqlist[v.POS-1+cur] = str(_alt)[cur]
            elif t=='INS': # is insertion
                newseqlist[v.POS-1] = str(_alt)
                # in case the REF is not a single base, take care of it
                for extra_i in range(_reflen-1):
                    curpos = v.POS+extra_i
                    if curpos not in newseqlist:
                        print("WARNING: {0}:{1} is already deleted! Check VCF format!".format(prefix, curpos))
                    else:
                        del newseqlist[curpos]
            else: # is deletion of size _d
                for i in range(abs(delta)):
                    curpos = v.POS+_reflen-2-i
                    if curpos not in newseqlist:
                        print("WARNING: {0}:{1} is already deleted! Check VCF format!".format(prefix, curpos))
                    else:
                        del newseqlist[curpos]

    vcf_writer.close()
    f_variant.close()

    f = open(output_fasta, 'w')
    newseq = make_seq_from_list(newseqlist, 0, len(refseq))
    f.write(">" + newid + "\n" + newseq + '\n')
    f.close()

    f = open(output_frag_fasta, 'w')
    i = 0
    j = 0  # init here, in the event that the entire sequence is 29903 "N"s, the frag.fasta file will be empty
    while i < len(newseqlist)-1 and newseqlist[i]=='N': i += 1
    while i < len(newseqlist)-1:
        # i is the first position that is not N
        j = i + 1  # j is now the second position that is not N in this segment
        # progress j until encountering the first N again, note some positions could be deleted, so ok to skip over them
        while j < len(newseqlist) and ((j not in newseqlist) or newseqlist[j]!='N'): j+=1
        f.write(">{0}_frag{1}\n{2}\n".format(newid, i+1, make_seq_from_list(newseqlist, i, j)))
        i = j + 1 # is now the second position that is N
        # progress i until encountering the first non-N again
        while i < len(newseqlist)-1 and ((i not in newseqlist) or newseqlist[i]=='N'): i+=1
    if j>i: f.write(">{0}_frag{1}\n{2}\n".format(newid, i+1, make_seq_from_list(newseqlist, i, j)))
    f.close()

    f = open(output_info, 'w')
    count_d = Counter(newseq.upper())
    f.write("total,num_A,num_T,num_C,num_C,num_N,num_sub,num_ins,num_del\n")
    f.write(str(len(newseq)) + ',' + \
            str(count_d['A']) + ',' + \
            str(count_d['T']) + ',' + \
            str(count_d['C']) + ',' + \
            str(count_d['G']) + ',' + \
            str(count_d['N']) + ',' + \
            str(tally_types['SUB']) + ',' + \
            str(tally_types['INS']) + ',' + \
            str(tally_types['DEL']) + '\n')
    f.close()


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("ref_fasta", help="Reference fasta (should be Wuhan ref)")
    parser.add_argument("prefix", help="Sample prefix")
    parser.add_argument("--input_depth", default=None, help="(optional) Input depth file, if not given, then <prefix>.bam.depth is expected.")
    parser.add_argument("--input_vcf", default=None, help="(optional) Input VCF file, if not given, then <prefix>.VCF is expected.")
    parser.add_argument("-c", "--min_coverage", type=int, default=4, help="Minimum base coverage to call a base (default: 4)")
    parser.add_argument("-f", "--min_alt_freq", type=float, default=0.5, help="Minimum variant frequency (default: 0.5)")
    parser.add_argument("-q", "--min_qual", type=int, default=100, help="Minimum QUAL cutoff (default: 100)")
    parser.add_argument("--vcf_type", required=True, choices=['pbaa', 'deepvariant', 'CLC', 'bcftools'], default=None, help="VCF format info")

    args = parser.parse_args()

    if args.min_alt_freq >= 1 or args.min_alt_freq <= 0:
        print("--min_alt_freq must be a fraction between (0,1]. Got {0} instead. Abort!".format(args.min_alt_freq))
        sys.exit(-1)

    if args.input_depth is None:
        depth_file = args.prefix + '.bam.depth'
    else:
        depth_file = args.input_depth

    if args.input_vcf is None:
        vcf_input = args.prefix + '.vcf'
    else:
        vcf_input = args.input_vcf

    if not os.path.exists(depth_file):
        print("Cannot find input file {0}. Abort!".format(depth_file))
        sys.exit(-1)

    if not os.path.exists(vcf_input):
        print("Cannot find input file {0}. Abort!".format(vcf_input))
        sys.exit(-1)

    prefix = args.prefix # ex: LC0003335
    newid = prefix + "_VCFConsensus"

    genVCFcons(args.ref_fasta, depth_file, vcf_input, args.prefix,
               newid=newid,
               min_coverage=args.min_coverage,
               min_alt_freq=args.min_alt_freq,
               min_qual=args.min_qual,
               vcf_type=args.vcf_type)
