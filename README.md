# CoSA
Coronavirus (SARS-Cov-2) sequencing analysis

Last Updated: 02.17.2020 (v6.1.0)

## Updates

02.17.2020    v6.1.0 release. Fixed `VCFCons.py` dealing with multiple ALT and proper INS if REF is more than 1nt.

02.15.2020    v6.0.0 release. Adding `consensusVariants.py` and `pbaa2vcf.py` for pbaa support.

02.14.2020    v5.2.0 release. `VCFCons.py` important bug fix on propagating deletions into consensus fasta.

02.11.2020    v4.0.0 release. `VCFCons.py` added.

## What is CoSA

CoSA is a set of Python and R scripts for analyzing SARS-CoV-2 sequences from PacBio HiFi/CCS data.  

See [wiki](https://github.com/Magdoll/CoSA/wiki) for details.


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
