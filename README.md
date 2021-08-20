# vcf2mst

Hamming Distance based Minimum Spanning Tree from Samples vcf using graptree
 - [Link to the Paper](https://www.medrxiv.org/content/10.1101/2021.05.25.21257370v1.article-metrics).

When using vcf2mst please use the following citation:

SARS-CoV-2 surveillance in Italy through phylogenomic inferences based on Hamming distances derived from functional annotations of SNPs, MNPs and InDels
Adriano Di Pasquale, Nicolas Radomski, Iolanda Mangone, Paolo Calistri, Alessio Lorusso, Cesare Camma
medRxiv 2021.05.25.21257370; doi: https://doi.org/10.1101/2021.05.25.21257370 

## Synopsis

```
usage 1: 
vcf2mst.pl samples_vcfcodes.csv > mst.nwk

usage 2: 
vcf2mst.pl list_of_vcfiles vcf > mst.nwk

```

See `examples` folder for `samples_vcfcodes.csv` and `list_of_vcfiles` file format

# Installation

**Prerequisites**

- Linux OS
- [Perl basic installation](https://www.perl.org/)
- [Graptree](https://github.com/achtman-lab/GrapeTree)

**Installation Step**

- clone/download vcf2mst on your system
- let scripts be executable, e.g., `chmod a+x *.pl`
- put perl scrypts in $PATH, e.g., `cp *.pl /usr/local/bin` 

## Docker based installation

- clone/download vcf2mst on your Linux based system
- `docker build -t vcf2mst . `

use it as in the following example

```
usage 1: 
docker run -u $UID -v /tmp:/tmp --rm  vcf2mst vcf2mst.pl samples_vcfcodes.csv > mst.nwk

usage 2: 
docker run -u $UID -v /tmp:/tmp --rm  vcf2mst vcf2mst.pl /tmp/list_of_vcfiles vcf > /tmp/mst.nwk 

```
