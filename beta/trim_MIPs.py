#!/usr/env/python3
from collections import defaultdict, namedtuple
from csv import DictWriter, DictReader
import pysam
import bisect
import pdb
import parasail


ProbeMatchInfo = namedtuple('ProbeMatchInfo', ['probename1', 'probename2', 'probename',
                                               'umi_len',
                                               'umi1', 'umi2',
                                               'shift1', 'shift2',
                                               'new_start', 'new_end',
                                               'size_diff',
                                               'samplename'])

HASHSIZE = 10
MAX_SW_PAD = 5
MIN_SW_SCORE = 20
SW_SCORE_MATRIX = parasail.matrix_create("ACGT", 2, -5)

def find_probe_matches_front(probeset, seq, ph_start, ph_end, max_shift=1):
    for shift in range(max_shift+1):
        q = seq[(ph_start+shift):(ph_end+shift)]
        if q in probeset:
            return shift, probeset[q]
    return None

def find_probe_matches_end(probeset, seq, ph_start, ph_end, max_shift=1):
    for shift in range(max_shift+1):
        q = seq[(-ph_end-shift):(-ph_start-shift)]
        if q in probeset:
            return shift, probeset[q]

def find_matching_probe_pair(m1, m2, readlen, readseq, proberegions):
    """

    :param m1: (shift, list of matching ext)
    :param m2: (shift, list of matching lig)
    :return: matched pairs if possible, otherwise whatever matches
    """
    # ex: (0, [('ML_COV_675f_22PPB_00080', 'ext', 'TTTGTCACGCACTCAAAGGGA')])
    shift1, list1 = m1
    shift2, list2 = m2
    match_candidates = []
    for name1, type1, fullseq1 in list1:
        if type1!='ext': continue # we want the 5' to always be ext, 3' always to be lig
        ext_start0, ext_end1, ext_seq = proberegions[name1]['ext']
        for name2, type2, fullseq2 in list2:
            if type2!='lig': continue
            lig_start0, lig_end1, lig_seq = proberegions[name2]['lig']
            expected_size = abs(lig_end1 - ext_start0)
            if name1==name2:
                return ((readlen-expected_size), name1, name2, shift1, shift2, fullseq1, fullseq2)
            match_candidates.append(((readlen-expected_size), name1, name2, shift1, shift2, fullseq1, fullseq2))
        # if we reach here, we have a 5' ext, but no matching 3' lig based on the hashing, so let's try to align it
        for name2, info in proberegions.items():
            lig_start0, lig_end1, lig_seq = info['lig']
            query_seq = readseq[-(len(lig_seq)+MAX_SW_PAD):]
            o1 = parasail.sg_qx_trace(query_seq, lig_seq, 3, 1, SW_SCORE_MATRIX)
            if o1.score >= MIN_SW_SCORE:
                expected_size = abs(lig_end1 - ext_start0)
                return ((readlen-expected_size), name1, name1, shift1, len(query_seq)-o1.end_query-1, ext_seq, lig_seq)

    # go through lig and find matching ext
    for name2, type2, fullseq2 in list2:
        if type2 != 'lig': continue
        lig_start0, lig_end1, lig_seq = proberegions[name2]['lig']
        for name1, info in proberegions.items():
            ext_start0, ext_end1, ext_seq = info['ext']
            query_seq = readseq[:(len(ext_seq)+MAX_SW_PAD)]
            o1 = parasail.sg_qx_trace(query_seq, ext_seq, 3, 1,
                                      SW_SCORE_MATRIX)
            if o1.score >= MIN_SW_SCORE:
                expected_size = abs(lig_end1 - ext_start0)
                return ((readlen - expected_size), name2, name2, len(query_seq)-o1.end_query-1, shift2, ext_seq, lig_seq)


    if len(match_candidates) == 0:
        return None
    else:
        # return match with minimal diff in expected size
        match_candidates.sort(key=lambda x: abs(x[0]))
        return match_candidates[0]


# def trim_MIPs(input_bam, output_prefix, probeset, proberegions, umi_len=5, fixed_trim_len=25, samplename='NA'):
#     """
#
#     :param input_bam: input CCS BAM
#     :param probeset: dict of (probe sequence) --> (probe name, "lig" or "ext")
#     :param proberegions: dict of (probe name) --> {'ext': (start0, end1), 'lig': (start0, end1)}
#     :param umi_len:  length of UMI flanking the probe sequences
#     :return:
#     """
#     ph_start = umi_len
#     ph_end = umi_len + HASHSIZE
#
#     reader = pysam.AlignmentFile(input_bam, 'rb', check_sq=False)
#     writer = pysam.AlignmentFile(output_prefix+'.bam', 'wb', template=reader)
#     f_info = open(output_prefix + '.info.txt', 'w')
#     writer_info = DictWriter(f_info,
#                              ['sample', 'read_id', 'read_len', 'diff_len', 'UMI_ext', 'UMI_lig',
#                               'extra_ext', 'extra_lig', 'probename_ext', 'probename_lig', 'strand'],
#                              delimiter='\t')
#     writer_info.writeheader()
#
#     for r in reader:
#         d = r.to_dict()
#         seqlen = len(r.seq)
#         #pdb.set_trace()
#         # try to determine strand based on probe matches
#         m1 = find_probe_matches_front(probeset, r.seq, ph_start, ph_end, max_shift=1)
#         m2 = find_probe_matches_end(probeset, r.seq, ph_start, ph_end, max_shift=1)
#         if m1 is not None or m2 is not None:
#             o = find_matching_probe_pair(m1 if m1 is not None else (None,[]),
#                                          m2 if m2 is not None else (None,[]),
#                                          seqlen, r.seq[ph_start:-ph_start], proberegions)
#             if o is not None:
#                 size_diff, probename1, probename2, shift1, shift2, fullseq1, fullseq2 = o
#                 size_diff  = size_diff - shift1 - shift2 - umi_len*2 # adjust size diff by removing UMI/extra len
#                 probename = probename1 + '--' + probename2
#             else:
#                 probename = 'Unassigned'
#                 probename1, probename2 = 'NA', 'NA'
#                 shift1 = 0
#                 shift2 = 0
#                 fullseq1 = 'N'*fixed_trim_len
#                 fullseq2 = fullseq1
#                 size_diff = 'NA'
#             # trim away the probes (for its full length), tag the barcode
#             new_start, new_end = ph_start+shift1+len(fullseq1), -(ph_start+shift2+len(fullseq2))
#             umi1 = r.seq[shift1:(shift1+umi_len)]
#             umi2 = r.seq[seqlen-(shift2+umi_len):seqlen-shift2]
#             d['seq'] = r.seq[new_start:new_end]
#             d['qual'] = d['qual'][new_start:new_end]
#             d['tags'].append('XM:Z:' + umi1 + umi2)
#             d['tags'].append('XC:Z:' + probename)
#             d['tags'].append('XA:Z:XM-XC')
#             x = pysam.AlignedSegment.from_dict(d, r.header)
#             writer.write(x)
#             # f_info.write("read_id\tUMI_ext\tUMI_lig\textra_ext\textra_lig\tprobename\n")
#             info = {'read_id': r.qname,
#                     'read_len': seqlen,
#                     'diff_len': size_diff,
#                     'UMI_ext': umi1,
#                     'UMI_lig': umi2,
#                     'extra_ext': r.seq[:shift1] if shift1>0 else 'NA',
#                     'extra_lig': r.seq[-shift2:] if shift2>0 else 'NA',
#                     'probename_ext': probename1,
#                     'probename_lig': probename2,
#                     'strand': '+',
#                     'sample': samplename}
#             writer_info.writerow(info)
#         else:
#             seq2 = str(Seq(r.seq).reverse_complement())
#             m1 = find_probe_matches_front(probeset, seq2, ph_start, ph_end, max_shift=1)
#             m2 = find_probe_matches_end(probeset, seq2, ph_start, ph_end, max_shift=1)
#             probename = 'Unassigned'
#             probename1, probename2 = 'NA', 'NA'
#             shift1 = 0
#             shift2 = 0
#             fullseq1 = 'N' * fixed_trim_len
#             fullseq2 = fullseq1
#             strand = '?'
#             size_diff = 'NA'
#             if m1 is not None or m2 is not None:
#                 strand = '-'
#                 o = find_matching_probe_pair(m1 if m1 is not None else (None, []),
#                                              m2 if m2 is not None else (None, []),
#                                              seqlen, seq2[ph_start:-ph_start], proberegions)
#                 if o is not None:
#                     size_diff, probename1, probename2, shift1, shift2, fullseq1, fullseq2 = o
#                     size_diff = size_diff - shift1 - shift2 - umi_len * 2 # adjust size diff by removing UMI/extra len
#                     probename = probename1 + '--' + probename2
#             # trim away the probes (for its full length), tag the barcode
#             new_start, new_end = ph_start + shift1 + len(fullseq1), -(ph_start + shift2 + len(fullseq2))
#             umi1 = seq2[shift1:(shift1 + umi_len)]
#             umi2 = seq2[seqlen - (shift2 + umi_len):seqlen - shift2]
#             d['seq'] = seq2[new_start:new_end]
#             d['qual'] = d['qual'][::-1]
#             d['qual'] = d['qual'][new_start:new_end]
#             d['tags'].append('XM:Z:' + umi1 + umi2)
#             d['tags'].append('XC:Z:' + probename)
#             d['tags'].append('XA:Z:XM-XC')
#             x = pysam.AlignedSegment.from_dict(d, r.header)
#             writer.write(x)
#             # f_info.write("read_id\tUMI_ext\tUMI_lig\textra_ext\textra_lig\tprobename\n")
#             info = {'read_id': r.qname,
#                     'read_len': seqlen,
#                     'diff_len': size_diff,
#                     'UMI_ext': umi1,
#                     'UMI_lig': umi2,
#                     'extra_ext': seq2[:shift1] if shift1 > 0 else 'NA',
#                     'extra_lig': seq2[-shift2:] if shift2 > 0 else 'NA',
#                     'probename_ext': probename1,
#                     'probename_lig': probename2,
#                     'strand': strand,
#                     'sample': samplename}
#             writer_info.writerow(info)
#     f_info.close()
#     writer.close()


def iter_cigar_string(cigar_string):
    i = 0
    while str.isnumeric(cigar_string[i]):
        i += 1
    num = cigar_string[:i]
    for s in cigar_string[i:]:
        if not str.isnumeric(s):
            yield int(num), s
            num = ''
        else:
            num += s

def iter_cigar_string_backwards(cigar_string):
    s = cigar_string[-1]
    i = len(cigar_string) - 2
    num = ''
    while i >= 0:
        if str.isnumeric(cigar_string[i]):
            num = cigar_string[i] + num
        else:
            yield int(num), s
            s = cigar_string[i]
            num = ''
        i -= 1
    yield int(num), s

def write_match_record(r, info, writer_bam, writer_info):
    d = r.to_dict()
    d['seq'] = r.seq[info.new_start:info.new_end]
    d['qual'] = d['qual'][info.new_start:info.new_end]
    d['tags'].append('XM:Z:' + info.umi1 + info.umi2)
    d['tags'].append('XC:Z:' + info.probename)
    d['tags'].append('XA:Z:XM-XC')
    d['cigar'] = '*'  # make it basically unmapped
    d['map_quality'] = '255'
    d['ref_name'] = '*'
    d['flag'] = '0'
    x = pysam.AlignedSegment.from_dict(d, r.header)
    writer_bam.write(x)
    # f_info.write("read_id\tUMI_ext\tUMI_lig\textra_ext\textra_lig\tprobename\n")
    info = {'read_id': r.qname,
            'read_len': len(d['seq']),
            'diff_len': info.size_diff,
            'UMI_ext': info.umi1,
            'UMI_lig': info.umi2,
            'extra_ext': r.seq[info.umi_len:info.shift1] if info.shift1 > info.umi_len else 'NA',
            'extra_lig': r.seq[-info.shift2:-info.umi_len] if info.shift2 > info.umi_len else 'NA',
            'probename_ext': info.probename1,
            'probename_lig': info.probename2,
            'strand': '-' if r.is_reverse else '+',
            'sample': info.samplename}
    writer_info.writerow(info)

SW_SEARCH_PAD = 10
def find_probe_match_by_mapping(r, seq2, proberegions, ext_regions_by_start, lig_regions_by_start):
    i = bisect.bisect_left(ext_regions_by_start, (r.reference_start - SW_SEARCH_PAD, '*'))
    lig_last_index = min(len(lig_regions_by_start)-1, bisect.bisect_right(lig_regions_by_start, (r.reference_end + SW_SEARCH_PAD, '*'), lo=i))
    while i < len(ext_regions_by_start):
        probename1 = ext_regions_by_start[i][1]
        ext_start0, ext_end1, ext_seq = proberegions[probename1]['ext']
        o1 = parasail.sg_qx_trace(seq2, ext_seq, 3, 1, SW_SCORE_MATRIX)
        if o1.score >= MIN_SW_SCORE:
            # now we need to know how long the UMI is
            num1, s1 = next(iter_cigar_string(str(o1.cigar.decode, 'utf-8')))
            if s1!='I':  # don't see the UMI! ugh!
                num1 = 0
            j = max(0, i - MAX_SW_PAD)
            while j <= lig_last_index:
                probename2 = lig_regions_by_start[j][1]
                lig_start0, lig_end1, lig_seq = proberegions[probename2]['lig']
                o2 = parasail.sg_qx_trace(seq2[-len(lig_seq) - MAX_SW_PAD:], lig_seq, 3, 1, SW_SCORE_MATRIX)
                if o2.score >= MIN_SW_SCORE:
                    #pdb.set_trace()
                    # it's a hit of (probename1, probename2)
                    num2, s2 = next(iter_cigar_string_backwards(str(o2.cigar.decode, 'utf-8')))
                    if s2!='I':
                        num2 = 0
                    expected_size = abs(lig_end1 - ext_start0)
                    size_diff = len(r.seq) - num1 - num2 - expected_size
                    return (size_diff, probename1, probename2,\
                            num1, num2, ext_seq, lig_seq)
                j += 1
        i += 1
    return None

def trim_MIPs_by_mapping(mapped_bam, output_prefix, proberegions, umi_len=5, fixed_trim_length=25, samplename='NA'):
    """
    Given the sequences already mapped to genome (mapped bam), identify the most likely probe combo and trim it
    If no probe pair is found, trim a fixed length based on <fixed_trim_length>
    """
    ext_regions_by_start = [] # sorted list of (start, probename) for ext
    lig_regions_by_start = []
    for probename, info in proberegions.items():
        ext_regions_by_start.append((info['ext'][0], probename))
        lig_regions_by_start.append((info['lig'][0], probename))

    # some quick sanity checking that proberegions is sorted
    for i in range(len(ext_regions_by_start)-1):
        assert ext_regions_by_start[i][0] < ext_regions_by_start[i+1][0]
        assert lig_regions_by_start[i][0] < lig_regions_by_start[i+1][0]

    reader = pysam.AlignmentFile(mapped_bam, 'rb', check_sq=False)
    writer = pysam.AlignmentFile(output_prefix+'.bam', 'wb', template=reader)
    f_info = open(output_prefix + '.info.txt', 'w')
    writer_info = DictWriter(f_info,
                             ['sample', 'read_id', 'read_len', 'diff_len', 'UMI_ext', 'UMI_lig',
                              'extra_ext', 'extra_lig', 'probename_ext', 'probename_lig', 'strand'],
                             delimiter='\t')
    writer_info.writeheader()

    already_assigned = set() # list of read_id that has been assigned
    unassigned_list = [] # use this to store unassigned, since sometimes a read can have multi-mappings

    for r in reader:
        if r.is_unmapped: continue
        if r.qname in already_assigned: continue

        m = find_probe_match_by_mapping(r, r.seq, proberegions, ext_regions_by_start, lig_regions_by_start)

        if m is None:
            # if we fail to find a match based on mapping, use the sequences themselves based on a hashed key
            flag_matched = False
            m1 = find_probe_matches_front(probeset, r.seq, ph_start=0, ph_end=HASHSIZE, max_shift=umi_len*2)
            m2 = find_probe_matches_end(probeset, r.seq, ph_start=0, ph_end=HASHSIZE, max_shift=umi_len*2)
            if m1 is not None or m2 is not None:
                o = find_matching_probe_pair(m1 if m1 is not None else (None, []),
                                             m2 if m2 is not None else (None, []),
                                             len(r.seq), r.seq, proberegions)
                if o is not None:
                    # we got a match!
                    size_diff, probename1, probename2, shift1, shift2, fullseq1, fullseq2 = o
                    #if probename1!=probename2: pdb.set_trace()
                    info = ProbeMatchInfo(probename1=probename1,
                                          probename2=probename2,
                                          probename=probename1 + '--' + probename2,
                                          umi_len=umi_len,
                                          umi1=r.seq[:shift1],
                                          umi2=r.seq[-shift2:],
                                          shift1=shift1,
                                          shift2=shift2,
                                          new_start=shift1 + len(fullseq1),
                                          new_end=len(r.seq)-shift2-len(fullseq2),
                                          size_diff=size_diff,
                                          samplename=samplename)
                    flag_matched = True
            if not flag_matched: # welp, still unassigned
                unassigned_list.append(r)
        else:
            flag_matched = True
            size_diff, probename1, probename2, shift1, shift2, ext_seq, lig_seq = m
            #if probename1 != probename2: pdb.set_trace()
            info = ProbeMatchInfo(probename1=probename1,
                                  probename2=probename2,
                                  probename=probename1 + '--' + probename2,
                                  umi_len=umi_len,
                                  umi1=r.seq[:min(umi_len, shift1)],
                                  umi2=r.seq[-min(umi_len, shift2):],
                                  shift1=shift1,
                                  shift2=shift2,
                                  new_start=shift1 + len(ext_seq),
                                  new_end=len(r.seq) - shift2 - len(lig_seq),
                                  size_diff=size_diff,
                                  samplename=samplename)

        if flag_matched:
            write_match_record(r, info, writer, writer_info)
            already_assigned.add(r.qname)

    # write all the unassigned here, checking they didn't exist in already_assigned
    for r in unassigned_list:
        if r.qname in already_assigned: continue
        info = ProbeMatchInfo(probename1='NA',
                              probename2='NA',
                              probename='Unassigned',
                              umi_len=umi_len,
                              umi1=r.seq[:umi_len],
                              umi2=r.seq[-umi_len:],
                              shift1=umi_len,
                              shift2=umi_len,
                              new_start=fixed_trim_length,
                              new_end=len(r.seq) - fixed_trim_length,
                              size_diff='NA',
                              samplename=samplename)
        write_match_record(r, info, writer, writer_info)
    f_info.close()
    writer.close()



def read_probes(probe_filename):
    """
    Read MIPs probe filename
    :return: probeset, a dict of hash --> (probe_id, ext/lig, sequence)
             proberegions, a dict of probe_id --> {'lig': (start0, end1), 'ext': (start0, end)}
    """
    # read the probe set
    probeset = defaultdict(lambda: [])
    proberegions = {} # probe name --> {'ext': (start0, end1), 'lig': (start0, end1)}
    reader = DictReader(open(probe_filename),delimiter='\t')
    for r in reader:
        s1, s2 = r['ext_probe_sequence'], r['lig_probe_sequence']
        probeset[s1[:HASHSIZE]].append((r['probe_id'], 'ext', s1))
        probeset[s2[-HASHSIZE:]].append((r['probe_id'], 'lig', s2))
        proberegions[r['probe_id']] = {'lig': (int(r['lig_probe_start(1start)'])-1, int(r['lig_probe_stop(1start)']), s2),
                                       'ext': (int(r['ext_probe_start(1start)'])-1, int(r['ext_probe_stop(1start)']), s1)}
    return probeset, proberegions

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("input_bam")
    parser.add_argument("output_prefix")
    parser.add_argument("probe_filename")
    parser.add_argument("--umi_len", default=5, type=int, help="UMI length (default: 5bp)")
    parser.add_argument("--fixed_trim_length", default=25, type=int, help="Default fixed trim length for those without a matching probe set (default: 25bp)")
    parser.add_argument("--samplename", default='NA', help="(optional) sample name")

    args = parser.parse_args()

    probeset, proberegions = read_probes(args.probe_filename)

    trim_MIPs_by_mapping(args.input_bam,
                         args.output_prefix,
                         proberegions,
                         umi_len=args.umi_len,
                         fixed_trim_length=args.fixed_trim_length,
                         samplename=args.samplename)

