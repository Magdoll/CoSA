#!/usr/bin/env python
import os, sys
from csv import DictReader, DictWriter
from collections import Counter, defaultdict
from collections import OrderedDict
from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np

prs_bcr = OrderedDict([('TGGGCC', 'IGG'), ('AATTCT', 'IGM'), ('AAGACA', 'IGK'), ('GGGAAC', 'IGL'), ('CGGGAA', 'IGA')])
bcr_mutation_dict = defaultdict(lambda: []) # subtype (IGG) --> mutated seqs {ex: TAGGCC, TGAGCC, etc}
for x in open('bcr_dict.csv'):
    mutseq, oriseq = x.strip().split(',')
    bcr_mutation_dict[prs_bcr[oriseq]].append(mutseq)

INFO_FIELDNAMES = ['id', 'strand', 'type', 'len', 'ilen', 'r2', 'umi', 'insert', 'primer', 'r1']
UMI_LEN = 12
LINKER = 'GTACGGG'  # Liz note: technically GTACGGGGG but I'm letting the last G slip

def identify_link_subclass(input_fastq, output_prefix, ambiguous_ok=False, find_prs_max=50):
    semi_good = []
    still_bad = []

    reader = SeqIO.parse(open(input_fastq),'fastq')
    for r in reader:
        s = str(r.seq)
        s2 = str(r.seq.reverse_complement())
        i, j = -1, -1
        for x,name in prs_bcr.items():
            i = s[:find_prs_max].find(x)
            if i > 0:
                print("precise found for ", x, " on + strand at pos", i)
                break
            j = s2[:find_prs_max].find(x)
            if j > 0: break

        if i < 0 and j < 0 and ambiguous_ok:
            for mutseq in bcr_mutation_dict[name]:
                i = s[:find_prs_max].find(mutseq)
                if i > 0:
                    print("ambiguous found for ", mutseq, " on + strand at pos", i)
                    break
        if i < 0 and j < 0 and ambiguous_ok:
            for mutseq in bcr_mutation_dict[name]:
                j = s2[:find_prs_max].find(mutseq)
                if j > 0:
                    print("ambiguous found for ", mutseq, " on + strand at pos", j)
                    break
        if i > 0: semi_good.append((r,len(r.seq)-i-6,name,'s1'))
        elif j > 0: semi_good.append((r,len(r.seq)-j-6,name,'s2'))
        else: still_bad.append((r,i))
        #if len(semi_good)>=1000: break

    seen = defaultdict(lambda: []) # (umi,type) --> list of CCS id
    seen_debug = defaultdict(lambda: Counter())  # (umi,type) --> (insert) --> count
    f_by_type = {}
    for rep_seq, rep_name in prs_bcr.items():
        f_by_type[rep_name] = open("{o}_{n}.determined.fq".format(o=output_prefix, n=rep_name), 'w')
    f2 = open(output_prefix + '.info.csv', 'w')
    f3 = open(output_prefix + '.cluster_info.csv', 'w')
    f3.write("tag\tcount\tmembers\n")
    fu = open(output_prefix + '.undetermined.fq', 'w')
    writer = DictWriter(f2, fieldnames=INFO_FIELDNAMES, delimiter='\t')
    writer.writeheader()

    #
    #  [r2 should mostly be blank] -- [12bp UMI] -- [insert] --- [primer] --- [r1, actually "C" region]
    #

    # for debugging
    linker_pos = []

    for p in semi_good:
        info = {'id': p[0].id, 'strand': p[-1], 'type': p[2], 'len': len(r.seq), 'ilen': 'NA', 'umi': 'NA', 'primer': 'NA', 'r2': 'NA', 'insert': 'NA', 'r1': 'NA'}
        if p[-1]=='s2':
            s = str(p[0].seq)
        else:
            s = str(p[0].seq.reverse_complement())
        i = s.find(LINKER)
        if i > 0: linker_pos.append(i)
        if UMI_LEN <= i < UMI_LEN*3:
            insert = s[i:p[1]]
            ilen = p[1]-i
            info['umi'] = s[i-12:i]
            info['r2'] = s[:i-12]
            info['primer'] = s[p[1]:p[1]+6]
            info['r1'] = s[p[1]+6:]
            info['insert'] = insert
            info['ilen'] = ilen
            tag = info['umi'] + '-' + info['type']
            writer.writerow(info)
            seen[tag].append(p[0].id)
            seen_debug[tag][insert] += 1
        else:
            SeqIO.write(p[0], fu, 'fastq')

    # now for each (umi,type), output the most common sequence
    umi_index = 0
    for umi_type in seen_debug:
        umi_index += 1
        umi, type = umi_type.split('-')
        major_seq, major_count = seen_debug[umi_type].most_common(1)[0]
        total_count = sum(seen_debug[umi_type].values())
        f_by_type[type].write("@pacbio.{0} UMI:{1}:{2} type:{3} mcount:{4} count:{5}\n".format(\
            umi_index, umi, 'G'*len(umi), type, major_count, total_count))
        f_by_type[type].write("{0}\n+\n{1}\n".format(major_seq, 'I'*len(major_seq)))

    for f in f_by_type.values(): f.close()
    f2.close()
    fu.close()

    for k,v in seen.items():
        f3.write("{0}\t{1}\t{2}\n".format(k, len(v), ",".join(v)))
    f3.close()

    linker_pos = np.array(linker_pos)
    print("DEBUG: # of linkers found", len(linker_pos))
    print("DEBUG: # of linkers found at 12bp:", sum(linker_pos==12))
    print("DEBUG: # of linkers found > 12bp:", sum(linker_pos>12))
    print("DEBUG: # of linkers found < 12bp:", sum(linker_pos<12))


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("fastq_filename", help="Input fastq, trimmed of read1/read2 primers")
    parser.add_argument("output_prefix", help="Output prefix")
    parser.add_argument("-n", "--ambiguous", default=False, action="store_true", help="Allow for ambiguous 6-bp subclass primer identification")
    parser.add_argument("-w", "--window", default=50, type=int, help="max window to search for 6-bp subclass primer (default: 50bp, use 100bp for IgGsub)")

    args = parser.parse_args()

    identify_link_subclass(args.fastq_filename, args.output_prefix, args.ambiguous, args.window)
