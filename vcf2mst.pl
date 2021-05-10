#!/usr/bin/env perl
#
#-------------------------------------
#   vcf2mst.pl 
#  
#   Hamming Distance based Minimum Spanning Tree from Samples vcf using graptree
#   usage 1: 
#   vcf2mst.pl samples_vcfcodes.csv > mst.nwk
#
#   usage 2: 
#   vcf2mst.pl list_of_vcfiles vcf > mst.nwk
#
#-------------------------------------

my ($f, $type)=@ARGV;

if($type eq 'vcf'){
    qx{
        vcflist2codes.pl $f > samples_vcfcodes.csv;
    }
    $f='samples_vcfcodes.csv';
};

qx{
    vcf2ham.pl $f > hdmatrix.tsv;
    graptree -p hdmatrix.tsv > mst.nwk 
}
