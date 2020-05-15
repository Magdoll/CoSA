# CoSA
Coronavirus (SARS-Cov-2) sequencing analysis

Last Updated: 05/15/2020 (v2.5.0)   

## Updates

2020.05.15 updated to v2.5.0. fixed `juliet_json_to_vcf.py` not handling multi-SNPs in a single codon.

2020.05.03 updated to v2.4.0. fixed again samtools fraction issue in `subsample_amplicons.py`

2020.04.16 updated to v2.1.0. fixed samtools fraction issue in `subsample_amplicons.py`

2020.04.15 updated to v2.0.0. Added scripts to support SARS-CoV-2 using PacBio HiFi/CCS data. See [wiki](https://github.com/Magdoll/CoSA/wiki) for details.

2020.04.01 updated to v1.1.0. Added fetching scripts for NCBI SRA.


## What is CoSA

CoSA is a set of Python and R scripts for analyzing SARS-CoV-2 sequences from:

* PacBio HiFi/CCS data.  See [wiki](https://github.com/Magdoll/CoSA/wiki) for details.
* From [GISAID](http://gisaid.org/) and [NCBI](https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/#nucleotide-sequences). 

This repository only hosts the code and the meta-analysis. 

To reproduce the GISAID results, you will need to be a registered GISAID user and download the fasta sequences and metadata CSV yourself.

## CoSA Report for GISAID 2020-03-31

[PDF](https://www.dropbox.com/s/is9flnbpn637ijx/gisaid_metadata_report.pdf?dl=0) version of the GISAID 2020-03-31 SARS-CoV-2 meta-analysis report.

The 2020-03-31 GISAID contains 2807 sequences, 2788 of which are human. 

![](https://github.com/Magdoll/CoSA/blob/master/latest_report/Rplot.top_country.png)

Below, the red arrows mark locations where there are consistent (more than 10+ sequences) stretches of "N"s, 
likely due to the same set of primer designs that fail to cover these regions?

![](https://github.com/Magdoll/CoSA/blob/master/latest_report/Screenshot%202020-04-01%2011.15.39.png)

## CoSA HowTo

* <a href="#req">GISAID data requirements</a>
* <a href="#install">How to install CoSA</a>
* <a href="#filter">How to clean up metadata and filter low-quality sequences</a>
* <a href="#report">How to produce the report PDF</a>
* <a href="#ncbi">Fetching NCBI sequences</a>

<a name="req"/>

### GISAID data requirements

After obtaining a [GISAID](http://gisaid.org/) account, you will need to 
download the fasta sequences and prepare a comma-separated metadata CSV file.

The fasta sequence IDs need to be in the format like the following (I'll add support for non-GISAID data later):
```
EPI_ISL_402119|hCoV-19/Wuhan/IVDC-HB-01/2019
EPI_ISL_402120|hCoV-19/Wuhan/IVDC-HB-04/2020
EPI_ISL_402121|hCoV-19/Wuhan/IVDC-HB-05/2019
EPI_ISL_402123|hCoV-19/Wuhan/IPBCAMS-WH-01/2019
```

The metadata CSV needs to be comma-separated, with at least the following columns:
```
Accession ID
Host
Specimen source
Sequencing technology
Location
```

<a name="install"/>

### How to install CoSA

The prerequisites are:
* Python 3.7
* [BioPython](https://biopython.org/)
* [R](https://www.r-project.org/), required only if you are generating the report figures

To install, clone the repo and install:

```
$ git clone https://github.com/Magdoll/CoSA.git
$ cd CoSA
$ python setup.py build
$ python setup.py install
```

<a name="filter"/>

### How to clean up metadata and filter low-quality sequences

**(1) Clean up metadata CSV**

Clean up the metadata CSV by running the script:
```
clean_up_metadata.py [metadata_csv]
```

which produces the cleaned output with the suffix `.modified.csv`.

**(2) Filter low quality sequences**

```
filter_gappedshort.py 
usage: [-m METADATA] [--min_length MIN_LENGTH] [--max_gaps MAX_GAPS]
       [--max_amb MAX_AMB]
       fasta_filename

positional arguments:
  fasta_filename        Input fasta filename

optional arguments:
  -h, --help            show this help message and exit
  -m METADATA, --metadata METADATA
                        Metadata CSV file (optional)
  --min_length MIN_LENGTH
                        Minimum sequence length (default: 28000)
  --max_gaps MAX_GAPS   Maximum stretches of 'N' gaps (default: 2)
  --max_amb MAX_AMB     Maximum number of ambiguous bases (default: 10)
```

For example, after running the two scripts on the fake dataset:

```
$ cd examples/
$ clean_up_metadata.py test.metadata.csv 
Output written to: test.metadata.modified.csv
$ filter_gappedshort.py test.fake.fasta -m test.metadata.modified.csv 
Output written to: test.fake.pass.fasta test.fake.fail.fasta test.fake.pass_fail.csv
``` 

The output csv file will contain additional columns that mark the sequences as FAIL/PASS.

<a name="report"/>
 
### How to produce the report PDF

Once you have the `.metadata.csv` or better, the `.pass_fail.csv`, run the R script to generate the report.

```
$ Rscript <path_to_CoSA>/R/summarize_metadata.R test.fake.pass_fail.csv
```

<a name="ncbi"/>

### Fetching NCBI sequences

You can use the following script to automatically fetch 
[NCBI SARS-CoV-2 GenBank records](https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/#nucleotide-sequences). 

The script will download the GenBank files and parse them into a fasta and 
metadata CSV that is the same format as the GISAID CSV I'm using in the rest of the tutorial.

```
fetch_NCBI.py examples/ncbi-2020-03-31.seqids.list ncbi_gbs/ 2020-03-31.ncbi
```

This will download individual GenBank filies into the directory `ncbi_gbs` 
(don't delete the directory, the script can auto-detect existing .gb files and skip re-downloading the next time)

The output will be collated into `2020-03-31.ncbi.fasta` and `2020-03-31.ncbi.metadata.csv`.
