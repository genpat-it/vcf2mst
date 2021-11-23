#!/usr/bin/env perl
#-------------------------------------
#   vcf2mst.pl 
#-------------------------------------
my $presentation=q{vcf2mst.pl         
Hamming Distance based Minimum Spanning Tree from Samples vcf using graptree

usage 1: 
vcf2mst.pl samples_vcfcodes.csv mst.nwk 

usage 2: 
vcf2mst.pl list_of_vcfiles mst.nwk  vcf

For additional info visit 
https://github.com/genpat-it/vcf2mst
};
my ($f, $out, $type)=@ARGV;

#-----------------------------------
# MAIN 
#-----------------------------------
init();

if( $type eq 'vcf'){
    $f=vcfListSnippy2Codes($f);
}elsif( $type eq 'gisaid'){
    $f=gisaidMetadata2Codes($f);
}

$f= vcf2ham($f);              #run("perl vcf2ham.pl $f > /tmp/hdmatrix.tsv");
run("$grapetreeCommand $f > $out");

print "DONE! newick file in $out \n";
#-----------------------------------



#-----------------------------------
# BASIC UTILS
#-----------------------------------
sub run{ my ($s) =@_;    
    print "$s";
    print "\nstart: " . qx{date};
    qx{$s};
    print "\nstop: " . qx{date};
}

#-----------------------------------
# vcfListSnippy2Codes
#-----------------------------------
sub vcfListSnippy2Codes{ my ($file) =@_;    
    #
    # Produce a tsv of sample\tvcfcode 
    # starting from vcf files in snippy format
    #
    my $ar=fileList2Array($file);

    my $out= "SAMPLECODE\tVCFCODE\n";
    my $err= '';
    foreach my $f (@{$ar}){
        if(! -e $f ){
            $err .= "# missing file: $f\n";
            next;
        }

        $cmp=$f;
        $cmp=~ s/^.*\///; # change the name of sample removing directories

        open(F, $f);
        while(<F>){
            chomp;
            if( $_ =~ /^#/){next;}

            my $r=$_;
            my @ar=split(/\t/,$r);
            my $type=''; my $code='';
            if( $r =~ /;TYPE=(\w+);/){
                $type=$1;
            }

            if( $type eq 'snp'){
                $code= "$type:@ar[3]@ar[1]@ar[4]";
            }elsif( $type eq 'mnp'){
                $code= "$type:@ar[3]@ar[1]@ar[4]";
            }elsif( $type =~ 'del' ){
                $pos=@ar[1] + 1;
                $pos2=$pos + length(@ar[3]) -2;
                $code= "$type:$pos-$pos2";
            }else{
                $code= "$type:@ar[1]-@ar[3]-@ar[4]";
            }
            $out .= "$cmp\t$code\n";
        }
    close(F);
    }
    #qx{         vcfListSnippy2Codes.pl $f > samples_vcfcodes.csv;     };
    
    my $out_file="/tmp/samples_vcfcodes.csv";
    open(F, ">$out_file" );
    print F $out;
    close(F);
    
    return $out_file;
}
#-----------------------------------

#-----------------------------
# fileList2Array
#-----------------------------
sub fileList2Array{ my ($f) =@_;
    my @ar;
    if( $f eq '' ){
        while(<>){
            chomp;
            push(@ar, $_);
        }
    }elsif( -e $f ){
        my $ris='';
        if( -d $f ){
            $ris=qx{find $f};
        }else{
            $ris=qx{cat $f};
        }
        @ar=split(/\n/,$ris);
	# use Data::Dumper;         print Dumper(\@ar);
    }else{
        die("file/dir not present: $f ")
    }

    return \@ar;
}

#-----------------------------------
# gisaidMetadata2Codes
#-----------------------------------

#
# usage 1: 
#   gisad2vcf.pl  samples_vcfcodes.csv 
#
#-------------------------------------
#


sub gisaidMetadata2Codes{ my ($file) =@_; 
    #
    # Produce a tsv of sample\tvcfcode 
    # starting from gisaid metadata file format
    # 
    # From
    #   hCoV-19/Tunisia/MHT_2/2020      (N_S202N,NS8_L84S)
    # to
    #   hCoV-19/Tunisia/MHT_2/2020      N_S202N
    #   hCoV-19/Tunisia/MHT_2/2020      NS8_L84S
    #
    my $out= "SAMPLECODE\tVCFCODE\n";
    open(F, $file);
    while(<F>){
        chomp;
        #print $_;
        if($_=~/^(\S+)[^\(]+\((\S+)\)/){
            $cmp=$1;
            $vcfstring=$2;
            #print "$cmp--- $vcfstring\n";
            @vcfs=split(/,/,$vcfstring);
            foreach $v (@vcfs) {
                $out .=  "$cmp\t$v\n";
            }
        }
    }
    close(F);
    
    my $out_file="/tmp/samples_vcfcodes.csv";
    open(F, ">$out_file" );
    print F $out;
    close(F);
    
    return $out_file;
}



#-------------------------------------
#   vcf2ham.pl 
#  
#   Hamming Distance matrix from Samples vcf codes
#   usage: 
#   vcf2mst.pl samples_vcfcodes.csv > hamming_distance_matrix.tsv
#
#   print a matrix compatible with grapetree input
#-------------------------------------
sub vcf2ham { my ($f)=@_;

    my $out_file="/tmp/hamming_distance_matrix.tsv";

    my $c={};
    my $codes={};
    my @a; 


    open(F, $f);
    while(<F>){
        chomp;
        if( $_=~ /SAMPLECODE/ ){next;}

        my ($SAMPLECODE,$VCFCODE) =   split(/\t/,$_);
        if( $VCFCODE eq ''){
            ($SAMPLECODE,$VCFCODE) =   split(/,/,$_);
        }
        
        if( ! exists($c->{$SAMPLECODE})){
                $c->{$SAMPLECODE}={};
        }

        $c->{$SAMPLECODE}->{$VCFCODE}=1;
        $codes->{$VCFCODE}++;
    }

    close(F);

    my $out="";
    #
    # out
    # 
    $out .= "#FILE\t";
    foreach $cod (sort(keys(%$codes))){    $out .= "$cod\t";}
    $out .= "\n";

    foreach $cmp (sort(keys(%$c))){
        $out .= "$cmp\t";
        foreach $cod (sort(keys(%$codes))){    
            my $val=( $c->{$cmp}->{$cod} )?1:2;
            $out .= "$val\t";
        }
        $out .= "\n";
    }

    open(F, ">$out_file" );
    print F $out;
    close(F);
    
    return $out_file;    
}

#--------------------------------------
#   INIT
#--------------------------------------

sub init{

    if( ! $f ){
        print $presentation;
        exit;          
    }

    if(! $type ){
        print "MODE: SAMPLECODE,VCFCODE \n";
    }elsif( $type =~ /^(vcf|gisaid)$/ ){
        print "MODE: $type \n";
    }else{
        print "ERROR: MODE NOT RECOGNIZED: $type\n";
        exit;
    }

    #$ENV{GRAPETREE_EXEC}='docker run  --mount type=bind,source=/tmp,destination=/tmp --rm quay.io/biocontainers/grapetree:2.1--pyh3252c3a_0 grapetree -p';

    if( ! $ENV{GRAPETREE_EXEC}){
        print "env GRAPETREE_EXEC is NOT set. grapetree will be used \n";
        $ENV{GRAPETREE_EXEC}='grapetree -p ';
    }else{
        print "env GRAPETREE_EXEC is set\n";
    }

    $grapetreeCommand=$ENV{GRAPETREE_EXEC};
    print "grapetreeCommand=$grapetreeCommand\n";

    if( ! $out ) {
        $out='/tmp/mst.nwk';
    };
}