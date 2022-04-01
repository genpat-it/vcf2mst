# vcf2mst

Hamming Distance based Minimum Spanning Tree from Samples vcf using graptree
 - [Link to the Paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-08112-0).

When using vcf2mst please use the following citation:


SARS-CoV-2 surveillance in Italy through phylogenomic inferences based on Hamming distances derived from functional annotations of SNPs, MNPs and InDels
Adriano Di Pasquale, Nicolas Radomski, Iolanda Mangone, Paolo Calistri, Alessio Lorusso, Cesare Camma
BMC Genomics 22, 782 (2021). https://doi.org/10.1186/s12864-021-08112-0


## Synopsis

```sh
#usage 1: 
vcf2mst.pl samples_vcfcodes.tsv mst.nwk code

#usage 2: 
vcf2mst.pl gisaid_metadata.tsv mst.nwk gisaid

#usage 3: 
vcf2mst.pl list_of_vcfiles mst.nwk  vcf

#usage 4: profile file only. return a matrix compatible with grapetree input
vcf2mst.pl samples_vcfcodes.tsv profile.tsv code   profile 
vcf2mst.pl list_of_vcfiles      profile.tsv vcf    profile 
vcf2mst.pl gisaid_metadata.tsv  profile.tsv gisaid profile 
```

See `examples` folder for `samples_vcfcodes.tsv` and `list_of_vcfiles` file format

See section **Usage** for further details  

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

- clone/download vcf2mst on your Linux/Unix based system
- `docker build -t vcf2mst . `


# Usage

## From vcfcodes csv file

```sh
# plain
vcf2mst.pl samples_vcfcodes.tsv mst.nwk vcf
# docker
docker run -u $UID -v /tmp:/tmp --rm  vcf2mst vcf2mst.pl /tmp/samples_vcfcodes.tsv  /tmp/mst.nwk
```

See `examples` folder for `samples_vcfcodes.tsv`  file format

The VCFCODE in samples_vcfcodes.tsv might be:

* a variant code from a vcf file in `snippy` format or
* a GISAID AA_SUBSTITUTIONS value from metadata.tsv (provided by GISAID) or 
* any other variant code produced by any `"variant caller"` of your own choice which maintain the condition `Same variant -> Same code` 

## From snippy format vcf files 

```sh
# plain
vcf2mst.pl list_of_vcfiles mst.nwk  vcf
# docker
docker run -u $UID -v /tmp:/tmp --rm  vcf2mst vcf2mst.pl /tmp/list_of_vcfiles  /tmp/mst.nwk vcf
```

See `examples/vcfiles` folder for  `list_of_vcfiles` file format

## From snippy format vcf director

```sh
# plain
vcf2mst.pl folder_with_vcfiles mst.nwk  vcf
# docker
docker run -u $UID -v /tmp:/tmp --rm  vcf2mst vcf2mst.pl /tmp/folder_with_vcfiles /tmp/mst.nwk vcf 
```

## From gisaid metadata.tsv 

```sh
# plain
vcf2mst.pl gisaid_metadata.tsv mst.nwk  gisaid
# docker
docker run -u $UID -v /tmp:/tmp --rm  vcf2mst vcf2mst.pl /tmp/gisaid_metadata.tsv /tmp/mst.nwk  gisaid
```
See `examples` folder for  `gisaid_metadata.tsv` file format (minimum information needed) and  `gisaid_full_metadata.tsv` (what is downloadable from gisaid) file format.


## GrapeTree command

In case the vcf2mst is installed locally and the docker version it is not used, there are different way to launch the required `grapetree` tool.
The environment variable GRAPETREE_EXEC can be used for different needs. 
In case there a local installation of grapetree is used, GRAPETREE_EXEC can be set like  `export GRAPETREE_EXEC=grapetree`  (this is the default)
In case there a docker installation of grapetree is used, GRAPETREE_EXEC can be set like  `export GRAPETREE_EXEC=docker run  --mount type=bind,source=/tmp,destination=/tmp --rm quay.io/biocontainers/grapetree:2.1--pyh3252c3a_0 grapetree -p`  

change the examples according to your specific installation and docker image


## Getting the "profile file"

the file 

