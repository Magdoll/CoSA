import vcf
import glob
import os, sys

file = sys.argv[1]

positions_B117 = [912, 3267, 5388, 5986, 6954, 14676, 15279, 16175, 21765, 21991, 23271, 23604, 23709, 24506, 24914, 27972, 28048, 28111, 28280, 28977]
positions_P1 = [733, 2838, 5648, 12778, 13860, 21621, 21638, 22810, 22811, 22812, 23010, 23011, 23012, 28167, 28262, 28512]

reader = vcf.VCFReader(open(file), 'rb')
for r in reader:
    if r.POS in positions_B117: print(file, r.POS, "B.1.1.7")
    elif r.POS in positions_P1: print(file, r.POS, "P1")
    
