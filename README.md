# CoSA
Coronavirus (SARS-Cov-2) sequencing analysis

Last Updated: 03/27/2020 (v1.0.0)   

## What is CoSA

CoSA is a set of Python and R scripts for analyzing SARS-CoV-2 sequences. 

For now, I'm only beginnig with [GISAID](http://gisaid.org/) data which hosts the largest number of submissions. I'm going to gradually look into adding other sources (ex: [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/)) and how to best combine multiple sources.

This repository only hosts the code and the meta-analysis. 
To reproduce the results, you will need to be a registered GISAID user and download the fasta sequences and metadata CSV yourself.

## CoSA Report for GISAID 2020-03-25

[PDF](https://www.dropbox.com/s/e1rqbk826uepdhz/gisaid_metadata_report.pdf?dl=0) version of the GISAID 2020-03-25 SARS-CoV-2 meta-analysis report.

The 2020-03-25 GISAID contains 1797 sequences, 1778 of which are human. All following analysis are done based on the human sequences.

![](https://github.com/Magdoll/CoSA/blob/master/latest_report/Screenshot%202020-03-27%2019.11.15.png?raw=true)

![](https://github.com/Magdoll/CoSA/blob/master/latest_report/Screenshot%202020-03-28%2005.59.45.png?raw=true)


Below, the red arrows mark locations where there are consistent (more than 10+ sequences) stretches of "N"s, likely due to the same set of primer designs that fail to cover these regions?

![](https://github.com/Magdoll/CoSA/blob/master/latest_report/Screenshot%202020-03-28%2006.04.33.png?raw=true)


## CoSA HowTo

* <a href="req">GISAID data requirements</a>
* <a href="install">How to install CoSA</a>
* <a href="filter">How to clean up metadata and filter low-quality sequences</a>
* <a href="report">How to produce the report PDF</a>

<a name="req"/>

### GISAID data requirements

After obtaining a [GISAID](http://gisaid.org/) account, you will need to 
download the fasta sequences and prepare a comma-separated metadata CSV file.

The fasta sequence IDs need to be in the format like the following:
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
