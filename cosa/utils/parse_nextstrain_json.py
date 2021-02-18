#!/usr/bin/env python
__author__ = "etseng@pacb.com"

import pdb
import json
from collections import namedtuple
from datetime import timedelta, datetime
from csv import DictReader, DictWriter

FIELDS_NUC = ['EPI', 'name', 'clade', 'mutation', 'date']
FIELDS_AA = ['EPI', 'name', 'clade', 'mutation', 'gene', 'date']

def convert_partial_year(number):
    year = int(number)
    d = timedelta(days=(number - year)*365)
    day_one = datetime(year,1,1)
    date = d + day_one
    return date.strftime("%Y-%m-%d")

def json_helper(node, writer_nuc, writer_aa):
    #pdb.set_trace()
    if 'children' not in node: # reached a tip
        # In [36]: x['branch_attrs']
        # Out[36]:
        # {'labels': {'aa': 'ORF1a: L2235I, N3833K'},
        #  'mutations': {'ORF1a': ['L2235I', 'N3833K'], 'nuc': ['C6968A', 'T11764A']}}
        epi = node['node_attrs']['gisaid_epi_isl']['value']
        name = node['name']
        clade = node['node_attrs']['clade_membership']['value']
        date = convert_partial_year(node['node_attrs']['num_date']['value'])
        if 'mutations' in node['branch_attrs']:
            for k,v in node['branch_attrs']['mutations'].items():
                if k == 'nuc':
                    for m in v:
                        writer_nuc.writerow({'EPI': epi, 'name': name, 'clade': clade, 'date': date, 'mutation': m})
                else:
                    for m in v:
                        writer_aa.writerow({'EPI': epi, 'name': name, 'clade': clade, 'date': date, 'mutation': m, 'gene': k})
    else:
        for c in node['children']:
            json_helper(c, writer_nuc, writer_aa)

def main(file):
    o = json.load(open(file))

    f_nuc = open(file + '.nuc_mutation.csv', 'w')
    f_aa = open(file + '.aa_mutation.csv', 'w')
    writer_nuc = DictWriter(f_nuc, FIELDS_NUC, delimiter=',')
    writer_aa = DictWriter(f_aa, FIELDS_AA, delimiter=',')
    writer_nuc.writeheader()
    writer_aa.writeheader()
    json_helper(o['tree'], writer_nuc, writer_aa)
    f_nuc.close()
    f_aa.close()

if __name__ == "__main__":
    file = 'ncov_global_2020-10-30_23-18.json'
    main(file)

