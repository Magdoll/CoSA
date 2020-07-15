import os, sys
from Bio import SeqIO

# @pacbio.{0} UMI:{1}:{2} type:{3} mcount:{4} count:{5}

def filter_by_umi_count(fastq_filename, cutoff):
    out_filename = fastq_filename[:fastq_filename.rfind('.')] + '.t' + str(cutoff) + '.fastq'
    f = open(out_filename, 'w')
    for r in SeqIO.parse(open(fastq_filename), 'fastq'):
        # ToDO: parse this using regex later
        raw = r.description.split('UMI:')
        umi = raw[1].split(':')[0]
        count = int(raw[1].split('count:')[1])
        if count >= cutoff:
            f.write("@{0} UMI:{1}:{2}\n{3}\n+\n{4}\n".format(r.id, umi, count, r.seq, 'I'*len(r.seq)))
    f.close()

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("fastq_filename")
    parser.add_argument("-c", "--cutoff", type=int, help="UMI count cutoff")

    args = parser.parse_args()

    filter_by_umi_count(args.fastq_filename, args.cutoff)
