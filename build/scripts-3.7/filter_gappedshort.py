#!/home/UNIXHOME/etseng/anacondaPy37/envs/anaCogentPy37/bin/python
"""
Filter step 1:
 - remove sequences that are too short
 - remove sequences with too many "N" stretches
"""

import os, sys, pdb, re
from csv import DictReader, DictWriter
from Bio import SeqIO
import pdb

seqid_rex = re.compile('(EPI_ISL_\d+)|(\S+)')

def trim_seq(seq):
    """
    Remove 'N's and non-ATCG from beginning and ends
    """
    good = ['A','T','C','G']
    seq = seq.upper()
    n = len(seq)
    i = 0
    while i < n and seq[i] not in good: i += 1
    j = len(seq)-1
    while j >= 0 and seq[j] not in good: j -= 1
    return seq[i:j+1]

def count_blocks(seq):
    """
    Yields (0-based start, 1-based end) of blocks of sequences not containing `N`s
    """
    n = len(seq)
    i = 0
    while seq[i]=='N': i += 1
    j = i + 1
    while j < n:
        if seq[j]=='N':
            yield (i,j)
            i = j + 1
            while i < n and seq[i]=='N': i += 1
            j = i + 1
        else:
            j += 1
    if i < n-1:
        yield (i,n)

def count_bad_base(seq):
    """
    Return the number of bases that are not A/T/C/G, excluding 'N's
    """
    #count = 0
    #for s in seq:
    #    if s not in ['A','T','C','G','N']:
    #        pdb.set_trace()
    #        print(s)
    return sum(s.upper() not in ['A','T','C','G','N'] for s in seq)

def filter_gappedshort(fasta_filename, min_len, max_gap, max_amb, csv_filename=None):

    output_prefix = fasta_filename[:fasta_filename.rfind('.')]
    f_pass = open(output_prefix + '.pass.fasta', 'w')
    f_fail = open(output_prefix + '.fail.fasta', 'w')
    f_csv  = open(output_prefix + '.pass_fail.csv', 'w')

    csv_fields = ['Accession ID', 'Length', 'GapCount', 'NCount', 'AmbBase', 'Pass1', 'Msg1']

    if csv_filename is None:
        writer = DictWriter(f_csv, csv_fields, delimiter=',')
    else:
        csv_reader = DictReader(open(csv_filename),delimiter=',')
        csv_info = {}
        for r in csv_reader:
            m = seqid_rex.match(r['Accession ID'])
            if m is None or m.group(1) is None:
                _id = r['Accession ID']
            else:
                _id = m.group(1)
            csv_info[_id] = r
        writer = DictWriter(f_csv, csv_reader.fieldnames + csv_fields[1:], delimiter=',')
    writer.writeheader()

    for r in SeqIO.parse(open(fasta_filename), 'fasta'):
        seq = trim_seq(str(r.seq))
        pass_flag, msg = True, 'NA'
        blocks = [(a, b) for (a, b) in count_blocks(seq)]
        num_gaps = len(blocks) - 1
        amb_base = count_bad_base(seq)

        if len(r.seq) < min_len:
            pass_flag, msg = False, "Too short"
        elif num_gaps > max_gap:
            pass_flag, msg = False, "Too many gaps"
        elif amb_base > max_amb:
            pass_flag, msg = False, "Too many ambiguous bases"
        if pass_flag:
            f_pass.write(">{0}\n{1}\n".format(r.id, seq))
        else:
            f_fail.write(">{0}\n{1}\n".format(r.id, seq))

        info = {'Accession ID': r.id,
                'Length': len(seq),
                'GapCount': num_gaps,
                'NCount': seq.count('N'),
                'AmbBase': amb_base,
                'Pass1': 'PASS' if pass_flag else 'FAIL',
                'Msg1': msg}
        if csv_filename is not None:
            m = seqid_rex.match(r.id)
            if m is None or m.group(1) is None:
                _id = r.id.split('|')[0]
            else:
                _id = m.group(1)
            info.update(csv_info[_id])

        stuff = ['"'+str(info[field])+'"' for field in writer.fieldnames]
        f_csv.write(",".join(stuff)+ '\n')

    f_pass.close()
    f_fail.close()
    f_csv.close()

    print("Output written to:", f_pass.name, f_fail.name, f_csv.name)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("")
    parser.add_argument("fasta_filename", help="Input fasta filename")
    parser.add_argument("-m", "--metadata", help="Metadata CSV file (optional)")
    parser.add_argument("--min_length", default=28000, type=int, help="Minimum sequence length (default: 28000)")
    parser.add_argument("--max_gaps", default=1, type=int, help="Maximum stretches of 'N' gaps (default: 1)")
    parser.add_argument("--max_amb", default=10, type=int, help="Maximum number of ambiguous bases (default: 10)")

    args = parser.parse_args()

    filter_gappedshort(args.fasta_filename, args.min_length, args.max_gaps, args.max_amb, args.metadata)