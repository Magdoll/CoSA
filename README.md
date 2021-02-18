# CoSA
Coronavirus (SARS-Cov-2) sequencing analysis

Last Updated: 02.17.2020 (v6.1.0)

[Stable version of CoSA under PacBio's GitHub](https://github.com/pacificbiosciences/CoSA)
[Developing version of CoSA under Magdoll's GitHub](https://github.com/Magdoll/CoSA)


## Updates

02.17.2020    v6.1.0 release. Fixed `VCFCons.py` dealing with multiple ALT and proper INS if REF is more than 1nt.

02.15.2020    v6.0.0 release. Adding `consensusVariants.py` and `pbaa2vcf.py` for pbaa support.

02.14.2020    v5.2.0 release. `VCFCons.py` important bug fix on propagating deletions into consensus fasta.

02.11.2020    v4.0.0 release. `VCFCons.py` added.

## What is CoSA

CoSA is a set of Python and R scripts for analyzing SARS-CoV-2 sequences from PacBio HiFi/CCS data.  

[Stable version Wiki](https://github.com/PacificBiosciences/CoSA/wiki)

[Developer version Wiki](https://github.com/Magdoll/CoSA/wiki)


<a name="install"/>

### How to install CoSA

To install, clone the repo and install:

```
$ git clone https://github.com/Magdoll/CoSA.git
$ cd CoSA
$ python setup.py build
$ python setup.py install
```
