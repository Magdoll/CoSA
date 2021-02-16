#!/usr/bin/env python
import re,pysam,os
import mappy as mp
import pandas as pd
import hashlib
from operator import itemgetter,xor
from collections import Counter
#from multiprocessing import Pool,cpu_count

#consensus sequence names
NAMEPATTERN=re.compile('sample-(?P<barcode>.*)_guide-(?P<guide>.*)_cluster-(?P<cluster>[0-9]+)_ReadCount-(?P<numreads>[0-9]+)')
DEFAULTSMAPNAMES=['Barcode','Bio Sample Name']
DEFAULTMINFRAC=0.01
DEFAULTPRESET='splice'
DEFAULTPREFIX='consensusVariants_summary'
DATEFORMAT='%Y-%m-%d %H:%M:%S'
ALLELETABLE='pbAA_consensus'
VARIANTTABLE='SampleVariants'

def main(parser):
    args = parser.parse_args()

    if args.noCSV and args.sqlite3 is None:
        raise ConsensusVariants_Error('No outputs defined with "--noCSV" and no sqlite3 db input')
        
    if xor(args.hifiSupport is None,args.read_info is None):
        raise ConsensusVariants_Error('Use both hifiSupport and read_info together')

    print('Loading Reference...')
    aligner = Aligner(args.reference,preset=None) #disabled presets at the moment
    
    sMap     = sampleMap(args.sampleMap)
    dtime    = args.datetime if args.datetime else getNow()
    alleles  = []
    variants = []
    print('Calling variants from consensus...')
    for consensusFa in args.consensusFastas:
        consensusType = parseHiLAAfastaName(consensusFa)
        for rec in pysam.FastxFile(consensusFa):
            if consensusType == 'failed' and nameDict.get('cluster_freq',0) < args.minFrac:
                continue #skip it

            aln = aligner(rec,skipFailed=(consensusType == 'failed')) 
            if aln is None:
                continue #skip failed consensus if they do not map
            nameDict  = parseName(rec)
            bioSample = sMap(nameDict['barcode'])
            if args.progress:
                print(f'Processing {bioSample}')
            #tableKey  = getKey(rec.name,bioSample,args.runName)
            tableKey  = getKey(rec.name,bioSample,args.runName,os.path.abspath(consensusFa))
            nameDict.update({'uuid'         :tableKey,
                             'runName'      :args.runName,
                             'bioSample'    :bioSample,
                             'source'       :os.path.abspath(consensusFa),
                             'clusterStatus':consensusType,
                             'chrom'        :aln.ctg,
                             'alnStart'     :aln.r_st,
                             'alnStop'      :aln.r_en,
                             'csString'     :aln.cs,
                             'length'       :len(rec.sequence),
                             'datetime'     :dtime})
            alleles.append(pd.Series(nameDict))
            variants.append(makeVarTable(aln,tableKey))

    #check for missing/failed results
    if len(sMap.remaining) and not args.ignoreMissing:
        for bc,sample in sMap.remaining:
            print(f'WARNINING: Failed sample {sample} / {bc}')
            tableKey = getKey(bc,sample,args.runName)
            alleles.append(pd.Series({'uuid'     :tableKey,
                                      'runName'  :args.runName,
                                      'bioSample':sample,
                                      'barcode'  :bc}))
    
    alleles_out = pd.DataFrame(alleles).set_index('uuid')
    variant_out = pd.concat(variants)

    if args.hifiSupport:
        print('Generating support alignments...')
        readSupport = supportEngine(args.read_info,
                                    args.hifiSupport,
                                    aligner)

        print('Appending reference calls...')
        variant_out = addRefCalls(variant_out,alleles_out)

        print('Tabulating support...')
        variant_out = readSupport(alleles_out,variant_out)

    if not args.noCSV:
        alleles_out.to_csv(f'{args.prefix}_alleles.csv')
        variant_out.to_csv(f'{args.prefix}_variants.csv')        
    if args.sqlite3:
        from sqlalchemy import create_engine
        engine = create_engine(f'sqlite:///{args.sqlite3}', echo=False)
        print(f'Importing to {args.sqlite3}')
        for pydf,sqldf in zip([alleles_out,variant_out],[ALLELETABLE,VARIANTTABLE]):
            pbaa = pydf.reset_index()\
                       .reindex(columns=list(TABLEMAP[sqldf].values()))
            pbaa.columns = list(TABLEMAP[sqldf].keys())
            pbaa.set_index('uuid').to_sql(sqldf, con=engine, if_exists='append')        

    return alleles_out,variant_out

class Aligner:
    def __init__(self,reference,preset=None):
        self.kwargs = {'fn_idx_in' : reference,
                       'best_n'    : 1}
        if preset:
            self.kwargs['preset']  = preset
        else:
            #self.kwargs['scoring'] = (2,5,5,4,56,0)
            self.kwargs['scoring'] = (1,2,2,1,32,0) # (A,B,o,e,O,E)
        self._aligner = mp.Aligner(**self.kwargs)
    
    def __call__(self,rec,skipFailed=True):
        try:
            return list(filter(lambda a:a.is_primary,self._aligner.map(seq=rec.sequence,cs=True)))[0]
        except IndexError:
            if skipFailed:
                return None
            else:
                raise ConsensusVariants_Error(f'Unable to align {rec.name}')

def getNow():
    from datetime import datetime
    now = datetime.now()
    return now.strftime(DATEFORMAT)

def makeVarTable(aln,key):
    variants = pd.DataFrame(parseCS(aln.cs,aln.r_st),columns=['POS','VAR'])
    variants['CHR'] = aln.ctg
    variants['uuid'] = key
    variants.set_index(['uuid','CHR','POS'],inplace=True)
    return variants


def addRefCalls(variants,alleles):
    #query to find ref-call consensus reads that span variants 
    qry = 'chrom==@ctg \
         & alnStart <= @pos <= alnStop \
         & uuid not in @vnts.index.get_level_values("uuid")'
    try:
        refCalls = pd.DataFrame([{'uuid': uuid,
                                  'CHR' : ctg,
                                  'POS' : pos,
                                  'VAR' : '.'}
                                 for (ctg,pos),vnts in variants.groupby(['CHR','POS'])
                                 for uuid in alleles.query(qry).index])\
                     .set_index((['uuid','CHR','POS']))    
        return pd.concat([variants,refCalls])
    except KeyError: #no records to add
        return variants

def safeFloat(v):
    try:
        return float(v)
    except ValueError:
        return v

def parseName(rec):
    resDict = NAMEPATTERN.search(rec.name).groupdict()
    splits  = list(map(lambda s:s.strip().split(),rec.comment.split(':')))
    flds    = map(itemgetter(-1),splits[:-1])
    vals    = [safeFloat(s[0]) for s in splits[1:-1]] + [' '.join(splits[-1])]
    resDict.update(dict(zip(flds,vals)))
    return resDict

def parseCS(csString,start=0,zeroIndex=True):
    ops = ':*-+~' 
    op  = None
    val = ''
    i   = start + int(zeroIndex)
    for s in csString:
        if s in ops:
            if op:
                incr,vnt = parseOp(i,op,val)
                if op=='+': #ins, place on prev position
                    i -= 1
                if vnt:
                    yield i,vnt
                i += incr
                val = ''
            op = s
        else:
            val += s
    incr,vnt = parseOp(i,op,val)
    if vnt:
        yield i,vnt

_splice = re.compile('[acgtn]{2}([0-9]+)[acgtn]{2}')
def parseOp(i,op,val):
    if op == ':': #match
        return int(val),None
    elif op == '-': #del
        return len(val),op+val
    elif op == '+': #ins
        return 1,op+val
    elif op == '*': #mismatch
        return 1,op+val
    elif op == '~': #splice/large del
        try:
            size = int(_splice.search(val).groups()[0])
        except Exception:
            raise ConsensusVariants_Error(f'Failed to parse variant: {op+val}')
        return size,f'{op}{size}del'
    else:
        raise ConsensusVariants_Error(f'Unknown variant type: {op}')

def parseHiLAAfastaName(fa):
    if fa.endswith('passed_cluster_sequences.fasta'):
        return 'passed'
    if fa.endswith('failed_cluster_sequences.fasta'):
        return 'failed'
    else:
        raise ConsensusVariants_Error(f'Input fasta {fa} not in HiLAA format')

def getKey(*args):
    return hashlib.md5(''.join(map(str,args)).encode()).hexdigest()

class supportEngine:
    FIELDNAME='supportReads'

    def __init__(self,readInfo,hifiReads,aligner):
        self.readInfo  = self._loadInfo(readInfo)
        self.alnGen    = self._align(hifiReads,aligner)
        self.vTable    = self._makeTable()
        #for now, only single sample runs
        assert self.vTable.Sample.nunique()==1, "can only do single samples"

    def __call__(self,alleles,variants):
        if variants.empty: #quick skip to avoid key error
            return variants.assign(**{self.FIELDNAME:None})

        sVar = pd.merge(alleles[['guide','cluster','numreads']],
                        variants,
                        left_index=True,
                        right_index=True)
        sVar.cluster = sVar.cluster.astype(int)

        bigTable = pd.merge(sVar.reset_index(),
                            self.vTable[['CHR','POS','VAR','guide','ClusterId']],
                            how='left',
                            left_on=['CHR','POS','guide','cluster'],
                            right_on=['CHR','POS','guide','ClusterId'],
                            suffixes=['_call','_reads'])

        counts   = bigTable.groupby(['uuid','CHR','POS'])\
                           .apply(self._countVars)\
                           .rename(self.FIELDNAME)

        return variants.join(counts)

    def _loadInfo(self,readInfo):
        return pd.read_csv(readInfo,sep=' ',names=READINFOCOLS)\
                 .set_index('readName')

    def _align(self,hifi,aligner):
        for rec in pysam.FastxFile(hifi):
            aln = aligner(rec)
            if aln is not None:
                yield rec.name,aln

    def _countVars(self,reads):
        numreads = int(reads.numreads.iloc[0])
        if reads.VAR_reads.isnull().all():
            return {'.':numreads}
        cnts     = Counter(reads.VAR_reads)
        subtot   = sum(cnts.values())
        if subtot < numreads:
            cnts['.'] += numreads - subtot
        return dict(cnts)

    def _makeTable(self):
        vTable = pd.concat([makeVarTable(aln,readName)
                            for readName,aln in self.alnGen])\
                   .reset_index(['CHR','POS'])

        #identify deletions >1 base
        #insert '*' if deleted positions are variants elsewhere
        #to prevent the reads from being counted as refcall over those pos
        multiDel     = vTable[vTable.VAR.str.startswith('-') & (vTable.VAR.str.len() > 2)]
        varPositions = set(map(tuple,vTable[['CHR','POS']].drop_duplicates().values))
        try:
            deletions    = pd.DataFrame([{'readName':row.name,
                                          'CHR'     :row.CHR,
                                          'POS'     :pos,
                                          'VAR'     :'*'}
                                          for i,row in multiDel.iterrows()
                                          for pos in range(row.POS+1,row.POS+len(row.VAR)-1)
                                           if (row.CHR,pos) in varPositions ])\
                                         .set_index('readName')
            return pd.concat([vTable,deletions]).join(self.readInfo)
        except KeyError:  #nothing to add
            return vTable.join(self.readInfo)

class sampleMap:
    def __init__(self,sampleMapFile=None):
        self.sMap      = {}
        self.remaining = []
        self.mapfile   = sampleMapFile
        #self.__call__  = self._getCallFunc(sampleMapFile)
        self.callFunc  = self._getCallFunc(sampleMapFile)

    def __call__(self,bc):
        return self.callFunc(bc)

    def _getCallFunc(self,mapFile):
        if mapFile is None:
            return (lambda bc: 'None')
        else:
            self.sMap      = dict(pd.read_csv(mapFile)[DEFAULTSMAPNAMES].values)
            self.remaining = list(self.sMap.items()) 
            return self._getSample
    _UNK='unknown_sample'
    def _getSample(self,bc):
        #sn = self.sMap[bc]
        sn = self.sMap.get(bc,self._UNK)
        if sn == self._UNK:
            print(f'WARNING:  Barcode {bc} not in sample map {self.mapfile}')
        try:
            self.remaining.pop(self.remaining.index((bc,sn)))
        except ValueError:
            pass
        return sn

class ConsensusVariants_Error(Exception):
    pass


TABLEMAP = {ALLELETABLE:      {'uuid'               :'uuid',
                               'runName'            :'runName',
                               'bioSample'          :'bioSample',
                               'barcode'            :'barcode',
                               'guide'              :'guide',
                               'cluster'            :'cluster',
                               'numreads'           :'numreads',
                               'uchime_score'       :'uchime_score',
                               'uchime_left_parent' :'uchime_left_parent',
                               'uchime_right_parent':'uchime_right_parent',
                               'cluster_freq'       :'cluster_freq',
                               'diversity'          :'diversity',
                               'avg_quality'        :'avg_quality',
                               'filters'            :'filters',
                               'source'             :'source',
                               'clusterStatus'      :'clusterStatus',
                               'alnStart'           :'alnStart',
                               'csString'           :'csString',
                               'length'             :'length',
                               'status'             :'clusterStatus',
                               'datetime'           :'datetime'},
            VARIANTTABLE:     {'uuid':'uuid',
                               'CHR' :'CHR',
                               'POS' :'POS',
                               'VAR' :'VAR'}
           }

READINFOCOLS = ['readName',
                'guide',
                'orient',
                'secondBestGuide',
                'Score',
                'ScoreParts',
                'Sample',
                'VarString',
                'ClusterId',
                'ClusterProb',
                'ClusterSize',
                'ChimeraScore']


if __name__ == '__main__':
    import argparse,sys

    parser = argparse.ArgumentParser(prog='consensusVariants.py', description='Write HiLAA consensus variants to table format')
    parser.add_argument('reference', metavar='reference', type=str,
                    help='Reference fasta or mmi')
    parser.add_argument('consensusFastas', metavar='consensusFastas', nargs='*', type=str,
                    help='Fasta file(s) of HiLAA consensus outputs')
    parser.add_argument('-r','--runName', dest='runName', type=str, default=None, required=True,
                    help=f'Sequencing run ID')
    parser.add_argument('-p','--prefix', dest='prefix', type=str, default=DEFAULTPREFIX, required=False,
                    help=f'Output prefix. Default {DEFAULTPREFIX}')
    parser.add_argument('-s','--sampleMap', dest='sampleMap', type=str, default=None,
                    help=f'CSV file mapping biosample name to barcode. Must have fields {DEFAULTSMAPNAMES}. Default None')
    parser.add_argument('-d','--datetime', dest='datetime', type=str, default=None,
                    help=f'Datetime of sequence run. Format "{DATEFORMAT.replace("%","%%")}". Default current time')
    parser.add_argument('-i','--ignoreMissing', dest='ignoreMissing', action='store_true', default=False,
                    help=f'Do not report empty rows for missing samples from sample map. Default False')
    parser.add_argument('-f','--minFrac', dest='minFrac', type=float, default=DEFAULTMINFRAC,
                    help=f'Ignore failed clusters below minFrac. Default {DEFAULTMINFRAC}')
    parser.add_argument('-P','--preset', dest='preset', choices=['splice','map-pb'], default=DEFAULTPRESET,
                    help=f'DISABLED. Alignment preset for mappy aligner. Choose "splice" for expected large deletions. Default {DEFAULTPRESET}')
    parser.add_argument('--hifiSupport', dest='hifiSupport', type=str, default=None,
                    help=f'Add per-variant read support depth. Path to hifi fastq used for clustering. Requires read_info option to be set. Default None')
    parser.add_argument('--read_info', dest='read_info', type=str, default=None,
                    help=f'Path to pbAA read_info table for the run. Required for hifi-support option. Default None')
    parser.add_argument('--sqlite3', dest='sqlite3', type=str, default=None,
                    help=f'Path to sqlite3 db to updload data. Will create or append to tables "{ALLELETABLE}" and "{VARIANTTABLE}". Default None')
    parser.add_argument('--noCSV', dest='noCSV', action='store_true', default=False,
                    help=f'Do not write CSV outputs. Default False (write to csv)')
    parser.add_argument('--progress', dest='progress', action='store_true', default=False,
                    help=f'Print sample names as they are processed. Default False')

    try:
        alleles,variants = main(parser)
    except ConsensusVariants_Error as e:
        print(f'\nERROR: {e}\n')
        sys.exit(1)



