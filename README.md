<a href="hhttps://programs.pacificbiosciences.com/l/1652/2021-03-12/42p2d7">![](https://github.com/Magdoll/images_public/blob/master/banners/test.png?raw=true)</a>

# CoSA
Coronavirus (SARS-Cov-2) sequencing analysis

Last Updated: 02.26.2020 (v8.5.0)

[Stable version of CoSA under PacBio's GitHub](https://github.com/pacificbiosciences/CoSA)

[Developing version of CoSA under Magdoll's GitHub](https://github.com/Magdoll/CoSA)

Join the [COVID19 Google Group](https://groups.google.com/g/smrt_covid19) to stay up-to-date on bioinfx recommendations and CoSA changes!
 

## Updates

02.26.2020    v8.5.0 release. cleaned up unnecesary directories.

02.23.2020    v8.4.0 release. `pbaa2vcf.py` format update fix.

02.22.2020    v8.3.0 release. `VCFCons.py ` now handles ill-formatted variants in CLC output.

02.21.2020    v8.2.0 release. `pbaa2vcf.py` now can handle edge cases of empty variants.

02.20.2020    v8.1.0 release. `VCFCons.py` can handle edge case of complete empty consensus and deals with unusual overlapping variant calls.

02.20.2020    v8.0.0 release. Added `.vcfcons.variants.csv` to `VCFCOns.py` output.

02.19.2020    v7.3.0 release. `VCFCons.py` accepts the new pbaa VCF type that follows standard AD info. Minor fixes to pbaa->VCF scripts.

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
$ git clone https://github.com/<ACCOUNT>/CoSA.git
$ cd CoSA
$ python setup.py build
$ python setup.py install
```
