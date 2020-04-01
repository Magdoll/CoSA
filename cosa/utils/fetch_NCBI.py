#!/usr/bin/env python
"""
Using NCBI eFetch to download and process GenBank into CSV files with the same metadata as the GISAID

Metadata:
File,Virus name,Accession ID,Type,Passage details/history,
Collection date,Location,Host,Additional location information,
Gender,Patient age,Patient status,Specimen source,
Additional host information,Outbreak,Last vaccinated,Treatment,
Sequencing technology,Assembly method,Coverage,Comment,
Originating lab,Address,Sample ID given by the sample provider,
Submitting lab,Address 1,Sample ID given by the submitting laboratory,
Authors,Submitter,Submission Date,Address 2,Continent,Country
"""

import os, sys
from csv import DictWriter
import urllib.request as rq
from Bio import GenBank
from Bio.SeqRecord import SeqRecord

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={id}&rettype=gb&retmode=text"
CSV_FIELDS = ['File', 'Virus name', 'Accession ID', 'Type', 'Passage details/history', 'Collection date',
              'Location', 'Host', 'Additional location information', 'Gender', 'Patient age', 'Patient status',
              'Specimen source', 'Additional host information', 'Outbreak', 'Last vaccinated',
              'Treatment', 'Sequencing technology', 'Assembly method', 'Coverage', 'Comment',
              'Originating lab', 'Address', 'Sample ID given by the sample provider', 'Submitting lab',
              'Address 1', 'Sample ID given by the submitting laboratory', 'Authors', 'Submitter',
              'Submission Date', 'Address 2', 'Continent', 'Country']

def process_host(x):
    """
    :return: host, gender, age
    """
    #Homo sapiens; male; age 65
    #Homo sapiens; female; age 49       Homo sapiens; male; age 41
    #Homo sapiens; female; age 52       Homo sapiens; male; age 61
    #Homo sapiens; male                 Homo sapiens; hospitalized patient
    raw = x.split('; ')
    gender = 'unknown'
    age = 'unknown'
    host = raw[0]
    if len(raw) == 3:
        if raw[1].lower() in ['male', 'female']: gender = raw[1].lower()
        if raw[2].startswith('age' ): age = raw[2].split('age ')[1]
    host = host.replace('"','')
    try:
        if host.lower().startswith('homo'): host = 'Human'
    except:
        pass
    return host, gender, age

def process_country(x):
    """
    :return: <continent/country>, <continent>, <country>
    """
    # China: Shenzhen
    # USA: WA
    raw = x.split(': ')
    if len(raw) == 1:
        raw = x.split(':')
    if len(raw) == 1:
        raw = x.split('/')
    if len(raw) == 1:
        country = raw[0]
        city = 'Unknown'
    else:
        country = raw[0]
        city = raw[1]

    country = country.replace('"', '')
    c = country.upper()
    if c == 'CHINA':
        continent = 'Asia'
        country = 'China'
    elif c == 'USA':
        continent = 'North America'
        country = 'USA'
    elif c.startswith("VIET"):
        continent = 'Asia'
        country = 'Vietnam'
    elif c in ['SPAIN', 'ITALY', 'BELGIUM', 'NETHERLANDS', 'SWEDEN']:
        continent = 'Europe'
    elif c in ['SOUTH KOREA', 'PHILIPPINES', 'THAILAND', 'TAIWAN' 'NEPAL', 'INDIA', 'JAPAN', 'PAKISTAN', 'IRAN']:
        continent = 'Asia'
    elif c in ['AUSTRALIA']:
        continent = 'Oceania'
    elif c in ['NIGERIA']:
        continent = 'Africa'
    elif c in ['BRAZIL', 'COLUMBIA']:
        continent = 'South/Central America'
    else:
        continent = 'Unknown'
    #print(c, country, continent, x)
    #input()
    return continent+'/'+country+'/'+city, continent, country


def process_genbank(gb_filename):
    """
    Process GenBank file name
    :param gb_filename: GenBank filename
    :return: dict of metadata, SeqRecord
    """
    info = {'File': None,
            'Virus name': None,
            'Accession ID': None,
            'Type': '',
            'Passage details/history': '',
            'Collection date': '',
            'Location': '',
            'Host': '',
            'Additional location information': '',
            'Gender': '',
            'Patient age' : '',
            'Patient status': '',
            'Specimen source': '',
            'Additional host information': '',
            'Outbreak': '',
            'Last vaccinated': '',
            'Treatment': '',
            'Sequencing technology': '',
            'Assembly method': '',
            'Coverage': '',
            'Comment': '',
            'Originating lab': '',
            'Address': '',
            'Sample ID given by the sample provider': '',
            'Submitting lab': '',
            'Address 1': '',
            'Sample ID given by the submitting laboratory': '',
            'Authors': '',
            'Submitter': '',
            'Submission Date': '',
            'Address 2': '',
            'Continent': '',
            'Country': ''}

    reader = GenBank.parse(open(gb_filename))
    r = next(reader)

    # In [55]: r.features[0]
    # Out[55]: Feature(key='source', location='1..29688')
    # In [56]: r.features[0].qualifiers
    # Out[56]:
    # [Qualifier(key='/organism=', value='"Severe acute respiratory syndrome coronavirus 2"'),
    #  Qualifier(key='/mol_type=', value='"genomic RNA"'),
    #  Qualifier(key='/isolate=', value='"SARS-CoV-2/WA-UW297/human/2020/USA"'),
    #  Qualifier(key='/host=', value='"Homo sapiens"'),
    #  Qualifier(key='/db_xref=', value='"taxon:2697049"'),
    #  Qualifier(key='/country=', value='"USA: WA"'),
    #  Qualifier(key='/collection_date=', value='"2020-03-15"')]

    info['File'] = gb_filename
    info['Accession ID'] = r.version
    info['Virus name'] = r.locus # later can rewrite if "/isolate=" is available
    info['Submission Date'] = r.date

    assert r.features[0].key=='source'
    q = dict((x.key,x.value) for x in r.features[0].qualifiers)
    if '/collection_date=' in q:
        info['Collection date'] = q['/collection_date=']
    if '/host=' in q:
        host, gender, age = process_host(q['/host='])
        info['Host'] = host
        info['Gender'] = gender
        info['Age'] = age
    if '/country=' in q:
        location, continent, country = process_country(q['/country='])
        info['Location'] = location
        info['Continent'] = continent
        info['Country'] = country
    if '/isolation_source=' in q:
        info['Specimen source'] = q['/isolation_source=']
    if '/isolate=' in q:
        info['Virus name'] = q['/isolate=']

    # In [94]: print(r.references[0])
    # REFERENCE   1  (bases 1 to 29838)
    #   AUTHORS   Chan,J.F.-W., Yuan,S., Kok,K.H., To,K.K.-W., Chu,H., Yang,J.,
    #             Xing,F., Liu,J., Yip,C.C.-Y., Poon,R.W.-S., Tsai,H.W., Lo,S.K.-F.,
    #             Chan,K.H., Poon,V.K.-M., Chan,W.M., Ip,J.D., Cai,J.P.,
    #             Cheng,V.C.-C., Chen,H., Hui,C.K.-M. and Yuen,K.Y.
    #   TITLE     A familial cluster of pneumonia associated with the 2019 novel
    #             coronavirus indicating person-to-person transmission: a study of a
    #             family cluster
    #   JOURNAL   Lancet (2020) In press
    #   REMARK    Publication Status: Available-Online prior to print
    info['Authors'] = r.references[0].authors
    info['Address'] = r.references[0].journal
    if len(r.references)>1:
        info['Submitter'] = r.references[1].authors
        info['Address 2'] = r.references[1].journal

    # In [113]: r.structured_comment
    # Out[113]:
    # OrderedDict([('Assembly-Data',
    #               OrderedDict([('Assembly Method', 'minimap2 v. 14 Jan 2020'),
    #                            ('Sequencing Technology', 'Nanopore')]))])
    try:
        if 'Assembly-Data' in r.structured_comment:
            x = r.structured_comment['Assembly-Data']
            if 'Sequencing Technology' in x:
                info['Sequencing technology'] = x['Sequencing Technology']
            if 'Assembly Method' in x:
                info['Assembly method'] = x['Assembly Method']
    except:
        pass

    # In [16]: r.taxonomy
    # Out[16]:
    # ['Viruses',
    #  'Riboviria',
    #  'Nidovirales',
    #  'Cornidovirineae',
    #  'Coronaviridae',
    #  'Orthocoronavirinae',
    #  'Betacoronavirus',
    #  'Sarbecovirus']
    info['Type'] = r.taxonomy[-2]

    _id = r.version
    if info['Virus name']!='':
        _id += '|' + info['Virus name']
    else:
        _id += '|' + r.definition
    return info, SeqRecord(r.sequence, id=_id)


def fetch_genbank_records(seqids, outdir, output_prefix):

    f_fasta = open(output_prefix+'.fasta', 'w')
    f_csv = open(output_prefix+'.metadata.csv', 'w')
    f_csv.write(",".join(CSV_FIELDS) + '\n')

    for seqid in seqids:
        gb_filename = os.path.join(outdir, seqid+'.gb')
        if not os.path.exists(gb_filename):
            print("Downloading {0}....".format(gb_filename), file=sys.stderr)
            url = EFETCH_URL.format(id=seqid)
            raw = rq.urlopen(url).read()
            with open(gb_filename, 'w') as f:
                f.write(raw.decode('utf-8'))
        print("Processing {0}....".format(gb_filename), file=sys.stderr)
        info, seqrec = process_genbank(gb_filename)

        f_fasta.write(">{0}\n{1}\n".format(seqrec.id, seqrec.seq))
        # stripping away nasty quotes before adding them back in a clean way
        for k in CSV_FIELDS:
            info[k] = info[k].replace('"' ,'')

        f_csv.write(",".join('"'+str(info[k])+'"' for k in CSV_FIELDS) + '\n')

    f_fasta.close()
    f_csv.close()


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("seq_list", help="Filelist of sequence IDs to download")
    parser.add_argument("out_dir", help="Output directory for downloaded .gb files")
    parser.add_argument("output_prefix", help="Output fasta/csv prefix")

    args = parser.parse_args()

    seqids = [line.strip() for line in open(args.seq_list)]
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    fetch_genbank_records(seqids, args.out_dir, args.output_prefix)






