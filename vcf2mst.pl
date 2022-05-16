#!/usr/bin/env perl
#-------------------------------------
#   vcf2mst.pl 
#-------------------------------------
my $presentation=q{vcf2mst.pl         
Hamming Distance based Minimum Spanning Tree from Samples vcf using graptree

#usage 1: 
vcf2mst.pl samples_vcfcodes.tsv  mst.nwk code

#usage 2: 
vcf2mst.pl gisaid_metadata.tsv   mst.nwk gisaid

#usage 3: 
vcf2mst.pl list_of_vcfiles       mst.nwk vcf

#usage 4: profile file only. return a matrix compatible with grapetree input
vcf2mst.pl samples_vcfcodes.tsv profile.tsv code   profile 
vcf2mst.pl list_of_vcfiles      profile.tsv vcf    profile 
vcf2mst.pl gisaid_metadata.tsv  profile.tsv gisaid profile 

For additional info visit 
https://github.com/genpat-it/vcf2mst

};

my ($f, $out, $type, $profile)=@ARGV;

#-----------------------------------
# MAIN 
#-----------------------------------
init();

# input
#
if( $type eq 'vcf'){
    $f=vcfListSnippy2Codes($f);
}elsif( $type eq 'gisaid'){
    $f=gisaidMetadata2Codes($f);
}elsif( $type eq 'code'){
}

# vcf2hammingdistance
#
$f= vcf2ham($f);             

if($profile){
    qx{cp $f $out};
    print "DONE! profile file in $out \n";
}else{
    #-----------------------------------
    # grapetree
    #-----------------------------------
    ###########profile2ids();
    run("$grapetreeCommand $f > $out");
    ###########ids2profile();    
    print "DONE! newick file in $out \n";
}
#-----------------------------------



#-----------------------------------
# BASIC UTILS
#-----------------------------------
sub run{ my ($s) =@_;    
    #
    #   runa a bash command taking start stop time
    #
    print "$s";
    print "\nstart: " . qx{date};
    qx{$s};
    print "\nstop: " . qx{date};
}

sub vcfListSnippy2Codes{ my ($file) =@_;    
    #-----------------------------------------------
    # Produce a tsv of sample\tvcfcode 
    # starting from vcf files in snippy format
    #-----------------------------------------------
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
            # split. fields from snippy            
            # CHROM	POS	TYPE	REF	ALT	EVIDENCE	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	LOCUS_TAG	GENE	PRODUCT
            my @ar=split(/\t/,$r); 
            #my $type=@ar[2]; 
            my ($chrom,$pos,$type,$ref,$alt)=@ar; 
            my $code='';
            if( $type eq 'snp'){
                $code= "$type:$ref$pos$alt";
            }elsif( $type eq 'mnp'){
                $code= "$type:$ref$pos$alt";
            }elsif( $type =~ 'del' ){
                $pos=$pos + 1;
                $pos2=$pos + length($ref) -2;
                $code= "$type:$pos-$pos2";
            }else{
                $code= "$type:$pos-$ref-$alt";
            }
            $out .= "$cmp\t$code\t$pos\n";
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

sub fileList2Array{ my ($f) =@_;
    #-----------------------------
    # fileList2Array
    #-----------------------------
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
    #-------------------------------------
    # Produce a tsv of sample\tvcfcode 
    # starting from gisaid metadata file format
    # 
    # From
    #   hCoV-19/Tunisia/MHT_2/2020      (N_S202N,NS8_L84S)
    # to
    #   hCoV-19/Tunisia/MHT_2/2020      N_S202N
    #   hCoV-19/Tunisia/MHT_2/2020      NS8_L84S
    #-------------------------------------
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
            my $pos;
            foreach $v (@vcfs) {
                #
                # pos calculation. codes examples
                # NSP3_I1413L,N_S33del,Spike_ins214EPE
                #   
                $pos='-1';
                if( $v=~/^(.*?_\w+\d+)/ ) {
                    $pos=$1;
                }
                $out .=  "$cmp\t$v\t$pos\n";
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
        if( $_=~ /TYPE/ ){next;}
        
        my ($SAMPLECODE,$VCFCODE,$POS) =   split(/\t/,$_);
        if( $VCFCODE eq ''){
            ($SAMPLECODE,$VCFCODE,$POS) =   split(/,/,$_);
        }
        
        if( ! exists($c->{$SAMPLECODE})){
            $c->{$SAMPLECODE}={};
        }

        ############### NEW ###############
        #
        # code based on position
        # 
        if( ! exists($codes->{$POS}) ){
            $max->{$POS}=0;
            $codes->{$POS}={};
        }
        if( ! exists($codes->{$POS}->{$VCFCODE}) ){
            my $allele_code = $max->{$POS} + 1;

            $codes->{$POS}->{$VCFCODE}=$allele_code;
            $max->{$POS}              =$allele_code;
        }
        
        $c->{$SAMPLECODE}->{$POS} = $codes->{$POS}->{$VCFCODE};
        
        ############### OLD ###############
        # $c->{$SAMPLECODE}->{$VCFCODE}=1;
        # $codes->{$VCFCODE}++;
        ###############
    }
    close(F);

    my $out="";
    #---------------------------
    # out
    #---------------------------

    #
    # Header
    #
    $out .= "#FILE\t";
    foreach $cod (sort(keys(%$codes))){    $out .= "$cod\t";}
    $out .= "\n";

    #
    # "allele" codes
    #
    foreach $cmp (sort(keys(%$c))){
        $out .= "$cmp\t";
        foreach $cod (sort(keys(%$codes))){
            my $allele=$c->{$cmp}->{$cod};
            my $val=( $allele )?
                    $allele + 1 :
                    1;
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
    }elsif( $type =~ /^(vcf|gisaid|code)$/ ){
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