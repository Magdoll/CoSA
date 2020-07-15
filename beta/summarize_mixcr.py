import os, sys
from csv import DictReader, DictWriter

#cloneId cloneCount      cloneFraction   clonalSequence  clonalSequenceQuality   allVHitsWithScore       allDHitsWithScore       all
#JHitsWithScore  allCHitsWithScore       allVAlignments  allDAlignments  allJAlignments  allCAlignments  nSeqFR1 minQualFR1      nSe
#qCDR1   minQualCDR1     nSeqFR2 minQualFR2      nSeqCDR2        minQualCDR2     nSeqFR3 minQualFR3      nSeqCDR3        minQualCDR3
#        nSeqFR4 minQualFR4      aaSeqFR1        aaSeqCDR1       aaSeqFR2        aaSeqCDR2       aaSeqFR3        aaSeqCDR3       aaS
#eqFR4   refPoints

OUT_FIELDNAMES = ['cloneId', 'cloneCount', 'cloneFraction', 'bestV', 'bestD', 'bestJ', 'bestC', 'clonalSequence',
                  'allVHitsWithScore', 'allDHitsWithScore', 'allJHitsWithScore', 'allCHitsWithScore']

def summarize_mixcr(mixcr_filename, output_filename):

    f = open(output_filename, 'w')
    writer = DictWriter(f, OUT_FIELDNAMES, delimiter='\t')
    writer.writeheader()

    for r in DictReader(open(mixcr_filename), delimiter='\t'):
        info = {'cloneId': r['cloneId'],
                'cloneCount': r['cloneCount'],
                'cloneFraction': r['cloneFraction'],
                'clonalSequence': r['clonalSequence'] if 'clonalSequence' in r else r['targetSequences'],
                'allVHitsWithScore': r['allVHitsWithScore'],
                'allDHitsWithScore': r['allDHitsWithScore'],
                'allJHitsWithScore': r['allJHitsWithScore'],
                'allCHitsWithScore': r['allCHitsWithScore']}
        # ex: IGHV4-59*00(3326),IGHV4-61*00(3115)
        info['bestV'] = r['allVHitsWithScore'].split(',')[0].split('*')[0]  # just keep the first one, ex: IGHV4-59
        info['bestD'] = r['allDHitsWithScore'].split(',')[0].split('*')[0]
        info['bestJ'] = r['allJHitsWithScore'].split(',')[0].split('*')[0]
        info['bestC'] = r['allCHitsWithScore'].split(',')[0].split('*')[0]
        writer.writerow(info)

    f.close()

if __name__ == "__main__":
    summarize_mixcr(sys.argv[1], sys.argv[1]+'.summary.txt')