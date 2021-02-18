#!/home/UNIXHOME/etseng/anacondaPy37/envs/anaCogentPy37/bin/python

from csv import DictReader

def assembly_method(x):
    x = x.upper()
    out = []
    if x.find('CLC GENOMICS')>=0:
        return 'CLC Genomics Workbench'
    elif x.find('ARTIC')>=0:
        return 'Artic pipeline'
    elif x.find('BWA')>=0:
        out.append('BWA')
    elif x.find('SAMTOOLS')>=0:
        out.append('Samtools')
    elif x.find('SPADES')>=0:
        out.append('SPAdes')
    elif x.find('MEGAHIT')>=0:
        out.append('MEGAHIT')
    elif x.find('GENEIOUS')>=0:
        out.append('Geneious ')
    elif x.find('MINIMAP')>=0:
        out.append('Minimap2')
    elif x.find('IVAR')>=0:
        out.append('ivar')
    elif x.find('NANOPOLISH')>=0:
        out.append('Nanopolish')
    elif x.find('BOWTIE2')>=0:
        out.append('Bowtie2')
    else:
        return "Unknown"

    out.sort()
    return ",".join(out)


def specimen_source(x):
    x = x.replace('-', '')
    x = x.split()
    x = " ".join(s[0].upper()+s[1:].lower() for s in x)
    if x.startswith('Naso') or x.startswith('Oraph') or x.startswith('Oroph') or x.startswith('Nose') or x.startswith('Throat') \
     or x.startswith('Nasa') or x.startswith('Phary') or x.startswith('Midnasal') or x.startswith('Mouth'):
        x = "NasoOropharyngeal Swab"
    elif x.startswith('Sputum'): x = 'Sputum'
    elif x.startswith('Alveolar') or x.startswith('Broncho'): x = 'Alveolar Lavage Fluid'
    elif x=='': x = 'Unknown'
    else: x = 'Other'
    return x
 
def location(x):
    raw = [t.strip() for t in x.split('/')]
    continent = raw[0]
    if len(raw)>=2: country = raw[1]
    else: country = 'Unknown'
    if country == 'China': continent = 'Asia'
    if continent.startswith('South America') or continent.startswith('Central Ameri'): continent = 'South/Central America'
    return continent, country

def tech(x):
    if x.find(';')>0 or x.find('+')>0 or x.find(',')>0: x = 'Multiple'
    elif x.startswith('Ion'): x = 'Ion Torrent'
    elif x.find('Illumina')>=0 or x.find('NextSeq')>=0 or x.find('Nova')>=0: x = 'Illumina'
    elif x.find('MGI')>=0: x = 'MGI'
    elif x.find('Sanger')>=0: x = 'Sanger'
    elif x.find('Nanopore')>=0 or x.find('nanopore')>=0 or x.find('ONT')>=0: x = 'ONT'
    elif x.find('Pac')>=0: x = 'PacBio'
    else: x = 'Unknown'
    return x


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("")
    parser.add_argument("csv_filename", help="Input CSV filename")

    args = parser.parse_args()

    input = args.csv_filename
    output = input[:input.rfind('.')] + '.modified.csv'

    f = open(output, 'w')
    reader = DictReader(open(input), delimiter=',')
    f.write(",".join(reader.fieldnames) + ',Continent,Country\n')
    for r in reader:
        stuff = []
        for x in reader.fieldnames:
            m = r[x]
            m = m.replace('"', '') # remove all double quotes
            if x=='Specimen source': m = specimen_source(r[x])
            elif x=='Sequencing technology': m = tech(r[x])
            elif x=='Assembly method': m = assembly_method(r[x])
            if m.startswith('"'): stuff.append(m)
            else: stuff.append('"'+m+'"')
        continent, country = location(r['Location'])
        f.write(",".join(stuff) + ',"' + continent + '","' + country + '"\n')

    f.close()
    print("Output written to:", f.name)