#import pdb
import os, sys, re
from collections import defaultdict, Counter
from Bio import SeqIO
import vcf


"""
Parser for `samtools mpileup`

http://www.htslib.org/doc/samtools-1.1.html
1. chr
2. 1-based position
3. ref base
4. coverage
5. readBase
6. base qualities
7. alignment qualities

readBase:
.  match to ref
,  match to ref on rev
> or <    ref skipping  (ex: like 37N)
ACGTN  mismatch on + strand
acgn   mismatch on - strand
+{number}{AGCTNagctn} - insertion of some {number}
-{number}{...} deletion of some {number}  # also means in next {number}, you will see a *
^ begin of read, followed by asci-33 for quality
$ end of read
"""


class MPileUpRecord(object):
    def __init__(self, chr, pos, ref, cov, readBase, baseQuals, alnQuals):
        """
        In addition to storing the 7 cols from mpileup,
        nalso stores
        counter: Counter of (key) -> (obs count in pileup)
        """
        self.chr = chr
        self.pos = pos
        self.ref = ref.upper() # let ref base always be upper case
        self.cov = cov
        self.nCov = None # this is the coverage of non-indel, non-skipped, which would be ACGTNacgtn
        self.nType = None # this is the number of non-indel, non-skipped bases accumulated at this record
        self.readBase = readBase
        self.baseQuals = baseQuals
        self.alnQuals = alnQuals

        self.counts = Counter()
        self.parse_readBase()

    def __str__(self):
        return """
        chr: {c}
        pos: {p} (1-based)
        ref: {r}
        cov: {v}
        nCov: {n}
        counts: {t}""".format(c=self.chr, p=self.pos+1, r=self.ref, v=self.cov, n=self.nCov, t=self.counts)

    def parse_readBase(self):
        """
        fill in self.counts
        """
        def not_indel_end_pos(i):
            return i >= len(self.readBase)-1 or self.readBase[i+1] not in ('+', '-', '$')

        rex = re.compile('(\d+)')
        def read_indel(start_index):
            m = rex.search(self.readBase, start_index)
            num = int(self.readBase[m.start():m.end()])
            return m.start(), m.end()+num

        sanity_counter = 0 # use this to track how many "reads" we've parsed to make sure parsing is correct
        # this number should agree with self.cov which is 4-th column in mpileup
        i = 0 # pointer for current location in string self.readBase
        while i < len(self.readBase):
            b = self.readBase[i]
            if b in '<>': # ignore skipped refs
                sanity_counter += 1
                i += 1
                continue
            elif b == '*': # deletion, just advance
                i += 1
                sanity_counter += 1
                continue
            elif b == '^': # start of read followed by ascii and either a comma or dot (ex: ^I.)
                i += 3
                sanity_counter += 1
                continue
            elif b == '$': # end of read, DO NOT advance counter
                i += 1
                continue
            elif b == '.': # could be followed by indels or $, careful don't double count
                self.counts[self.ref] += 1
                sanity_counter += 1
                i += 1
            elif b == ',': # # could be followed by indels or $, careful don't double count
                self.counts[self.ref.lower()] += 1
                sanity_counter += 1
                i += 1
            elif b in 'ATCGNatcgn':
                self.counts[b] += 1
                sanity_counter += 1
                i += 1
            elif b == '-': # DO NOT ADVANCE the sanity counter! otherwise double counting
                start, end = read_indel(i+1)
                self.counts["-"+self.readBase[start:end]] += 1
                i = end
            elif b == '+': # insertion should be +{number}{bases}
                start, end = read_indel(i+1)
                self.counts["+"+self.readBase[start:end]] += 1
                i = end
            else:
                raise Exception("Unknown {0} in readBase!".format(b))

        assert self.cov == sanity_counter or (self.readBase=='*' and self.cov==0)
        # set nCov which is cov provided by non-indel non-skipped bases
        self.nCov = 0
        self.nType = 0
        for x in 'ATCGNatcgn':
            self.nCov += self.counts[x]
            if self.counts[x] > 0: self.nType += 1


class MPileUpReader(object):
    def __init__(self, filename):
        self.filename = filename
        self.f = open(filename)

    def __iter__(self):
        return self

    def __next__(self):
        cur = self.f.tell()
        line = self.f.readline()
        if self.f.tell() == cur:
            raise StopIteration
        return self.parseLine(line)

    def parseLine(self, line):
        raw = line.strip().split('\t')
        if (len(raw)==7 or len(raw)==15):
            cov = int(raw[3])
            #if cov > 0:
            return MPileUpRecord(chr=raw[0], \
                                 pos=int(raw[1])-1, \
                                 ref=raw[2],
                                 cov=int(raw[3]),
                                 readBase=raw[4],
                                 baseQuals=raw[5],
                                 alnQuals=raw[6])
        elif len(raw)==4:
            # only way to have only 4 columns is because after --min-BQ filtering there are no bases
            # ex:
            # fake    8728    T       3       .$.$.   ;q:     ]]]
            # fake    8729    T       0
            return MPileUpRecord(chr=raw[0], \
                                 pos=int(raw[1])-1, \
                                 ref=raw[2],
                                 cov=0,
                                 readBase='',
                                 baseQuals='',
                                 alnQuals='')
        else:
            raise Exception("Expected to have 7 cols in mpileup record \
            but saw only {0}, abort! Line was: {1}".format(len(raw), line))


def make_seq_from_list(seqlist, start0, end1):
    seq = ''
    for i in range(start0, end1):
        if i in seqlist:
            seq += seqlist[i]
    return seq

def genVCFcons(ref_fasta, mpileup, vcf_input, prefix, newid, min_coverage=4, min_alt_freq=0.5, min_qual=100):
    """

    :param ref_fasta: should be the Wuhan reference
    :param mpileup: <sample>.bam.mpileup of aligning to Wuhan ref
    :param vcf_input: VCF of where the variants are
    :param prefix: output prefix
    :param min_coverage: below this coverage bases will be 'N'
    :param min_alt_freq: below this ALT frequency bases will use the reference instead
    :return:
    """
    #mpileup = prefix + '.bam.mpileup'
    #vcf_input = prefix + '.vcf'
    output = prefix + '.vcfcons.fasta'
    output2 = prefix + '.vcfcons.frag.fasta'

    ref = next(SeqIO.parse(open(ref_fasta),'fasta'))
    refseq = str(ref.seq)
    refseq = list(refseq)

    reader = MPileUpReader(mpileup)
    info_per_pos = {} # 0-based POS --> mpileup record
    for r in reader: info_per_pos[r.pos] = r

    newseqlist = dict(zip(range(len(refseq)), list(refseq)))
    for pos0 in range(len(refseq)):
        if pos0 not in info_per_pos or info_per_pos[pos0].cov < min_coverage: newseqlist[pos0] = 'N'
    # make sure begin/ends are "N"s
    for pos0 in range(min(info_per_pos)): newseqlist[pos0] = 'N'
    for pos0 in range(max(info_per_pos),len(refseq)): newseqlist[pos0] = 'N'

    # now add in the variants
    vcf_reader = vcf.Reader(open(vcf_input))
    vcf_writer = vcf.Writer(open(prefix+'.vcfcons.vcf', 'w'), vcf_reader)
    for v in vcf_reader:
        if len(v.ALT)>1:
            print("WARNING: more than 1 alt type for {0}! Using just the first alt for now.".format(prefix))
            #sys.exit(-1)

        _ref, _alt = str(v.REF), str(v.ALT[0]) # we'll ignore multi variants for now
        delta = len(_alt)-len(_ref)
        if delta==0: t = 'SUB'
        elif delta>0: t = 'INS'
        else: t = 'DEL'

        mrec = info_per_pos[v.POS-1]
        if t=='SUB':
            alt_count = mrec.counts[_alt] + mrec.counts[_alt.lower()]
        elif t=='INS':
            tmp = '+' + str(delta)
            alt_count = mrec.counts[tmp+_alt[1:]] + mrec.counts[tmp+_alt[1:].lower()]
        else:
            tmp = str(delta)  # delta is already -<something>, don't need to add '-' sign
            alt_count = mrec.counts[tmp+_ref[1:]] + mrec.counts[tmp+_ref[1:].lower()]

        alt_freq = alt_count * 1. / mrec.cov
        if alt_freq < min_alt_freq:
            print("WARNING: For {0}: Ignore variant {1}:{2}->{3} because alt freq is {4}.".format(prefix, v.POS, _ref, _alt, alt_freq))
        elif v.QUAL < min_qual:
            print("WARNING: For {0}: Ignore variant {1}:{2}->{3} because qual is {4}.".format(prefix, v.POS, _ref, _alt, v.QUAL))
        else:
            vcf_writer.write_record(v)
            if t=='SUB':
                newseqlist[v.POS-1] = str(v.ALT[0])
            elif t=='INS': # is insertion
                newseqlist[v.POS-1] = str(v.ALT[0])
            else: # is deletion of size _d
                for i in range(delta): del newseqlist[v.POS+i]

    vcf_writer.close()
    f = open(output, 'w')
    newseq = make_seq_from_list(newseqlist, 0, len(refseq))
    f.write(">" + newid + "\n" + newseq + '\n')
    f.close()

    f = open(output2, 'w')
    i = 0
    while newseqlist[i]=='N': i += 1
    while i < len(newseqlist)-1:
        j = i + 1
        while j < len(newseqlist) and newseqlist[j]!='N': j+=1
        f.write(">{0}_frag{1}\n{2}\n".format(newid, i+1, make_seq_from_list(newseqlist, i, j)))
        i = j + 1
        while i < len(newseqlist)-1 and newseqlist[i]=='N': i+=1
    if j>i: f.write(">{0}_frag{1}\n{2}\n".format(newid, i+1, make_seq_from_list(newseqlist, i, j)))
    f.close()


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("ref_fasta", help="Reference fasta (should be Wuhan ref)")
    parser.add_argument("prefix", help="Sample prefix")
    parser.add_argument("-s", "--seq_rename_filename", default=None, help="Sequence ID rename file name. Optional.")
    parser.add_argument("-c", "--min_coverage", type=int, default=4, help="Minimum base coverage to call a base (default: 4)")
    parser.add_argument("-f", "--min_alt_freq", type=float, default=0.5)
    parser.add_argument("-q", "--min_qual", type=int, default=100, help="Minimum QUAL cutoff (default: 100)")

    args = parser.parse_args()

    if args.min_alt_freq >= 1 or args.min_alt_freq <= 0:
        print("--min_alt_freq must be a fraction between (0,1]. Got {0} instead. Abort!".format(args.min_alt_freq))
        sys.exit(-1)

    mpileup = args.prefix + '.bam.mpileup'
    vcf_input = args.prefix + '.vcf'

    if not os.path.exists(mpileup):
        print("Cannot find input file {0}. Abort!".format(mpileup))
        sys.exit(-1)

    if not os.path.exists(vcf_input):
        print("Cannot find input file {0}. Abort!".format(vcf_input))
        sys.exit(-1)

    prefix = sys.argv[1] # ex: LC0003335


    #read CDC's renaming file
    newid = None #args.prefix+'_VCFconsensus'


    if args.seq_rename_filename is not None:
        for line in open(args.seq_rename_filename):
            # ex: hCoV-19/USA/IA-CDC-LC0005111/2021
            tmp = line.strip().split('/')[-2].split('-')[-1]
            if tmp == args.prefix:
                newid = line.strip()
    if newid is None:
        newid = args.prefix+'_VCFconsensus'

    genVCFcons(args.ref_fasta, mpileup, vcf_input, args.prefix, newid, min_coverage=args.min_coverage, min_alt_freq=args.min_alt_freq, min_qual=args.min_qual)