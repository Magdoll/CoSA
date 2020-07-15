import os, sys
from csv import DictWriter

SUMMARY_HEAD = "V-(D)-J rearrangement summary for query sequence"
FIELDS = ['pbid', 'length', 'Vmatch', 'Dmatch', 'Jmatch', 'VmatchTop', 'DmatchTop', 'JmatchTop', 'chaintype', 'stopcodon', 'VJframe', 'productive', 'strand']

def process_igblast(ig_filename):

    f = open(ig_filename + '.summary.csv', 'w')
    writer = DictWriter(f, FIELDS, delimiter='\t')
    writer.writeheader()
    flag_vdjsummary = False

    for line in open(ig_filename):
        if line.startswith('Query= '):
            pbid = line.strip().split('Query= ')[1]
        elif line.startswith('Length='):
            pblen = line.strip().split('Length=')[1]
        elif line.startswith(SUMMARY_HEAD):
            flag_vdjsummary = True
        elif flag_vdjsummary:
            raw = line.strip().split()
            if len(raw) == 8:
                info = {'pbid': pbid,
                        'length': pblen,
                        'Vmatch': raw[0],
                        'VmatchTop': raw[0].split(',')[0],
                        'Dmatch': raw[1],
                        'DmatchTop': raw[1].split(',')[0],
                        'Jmatch': raw[2],
                        'JmatchTop': raw[2].split(',')[0],
                        'chaintype': raw[3],
                        'stopcodon': raw[4],
                        'VJframe': raw[5],
                        'productive': raw[6],
                        'strand': raw[7]}
            elif len(raw) == 7:
                info = {'pbid': pbid,
                        'length': pblen,
                        'Vmatch': raw[0],
                        'VmatchTop': raw[0].split(',')[0],
                        'Dmatch': 'NA',
                        'DmatchTop': 'NA',
                        'Jmatch': raw[1],
                        'JmatchTop': raw[1].split(',')[0],
                        'chaintype': raw[2],
                        'stopcodon': raw[3],
                        'VJframe': raw[4],
                        'productive': raw[5],
                        'strand': raw[6]}
            else:
                print("{0} has incorrect number of entries!".format(len(raw)))
                sys.exit(-1)
            writer.writerow(info)
            flag_vdjsummary = False
            pbid, pblen = None, None
    f.close()

process_igblast(sys.argv[1])





