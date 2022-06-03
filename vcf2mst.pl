#!/usr/bin/env perl
#-------------------------------------
#   vcf2mst.pl 
#-------------------------------------
my $presentation=q{vcf2mst.pl         
Hamming Distance based Minimum Spanning Tree from Samples vcf using graptree.

 - [Link to the Paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-08112-0).

When using vcf2mst please use the following citation:

SARS-CoV-2 surveillance in Italy through phylogenomic inferences based on Hamming distances derived from functional annotations of SNPs, MNPs and InDels
Adriano Di Pasquale, Nicolas Radomski, Iolanda Mangone, Paolo Calistri, Alessio Lorusso, Cesare Camma
BMC Genomics 22, 782 (2021). https://doi.org/10.1186/s12864-021-08112-0

## Synopsis

```sh
# Basic usage
vcf2mst.pl input_tsv_file out_file type_of_input [options]

# type_of_input=[vcf,gisaid,algn2pheno,nextclade,tsv,code]

# Examples

# type of inputs
vcf2mst.pl list_of_vcfiles           mst.nwk vcf
vcf2mst.pl gisaid_metadata.tsv       mst.nwk gisaid
vcf2mst.pl algn2pheno_metadata.tsv   mst.nwk algn2pheno
vcf2mst.pl nextclade_metadata.tsv    mst.nwk nextclade
vcf2mst.pl samples_vcfcodes.tsv      mst.nwk code

# profile file only. return a matrix compatible with grapetree input
vcf2mst.pl list_of_vcfiles      profile.tsv vcf    -out profile 

# usage 5: filter positions
vcf2mst.pl list_of_vcfiles      profile.tsv vcf    -minNT 10 -maxNT 10000
```

For additional info visit 
https://github.com/genpat-it/vcf2mst

};
my ($f, $out, $type, $options)=@ARGV;
my %opt;
#-----------------------------------
# MAIN 
#-----------------------------------
init();

#
# input 2 variant codes
#
if( $type eq 'vcf'){                $f=vcfListSnippy2Codes($f);
}elsif( $type =~ /(nextclade|algn2pheno|gisaid)/ ){     
                                    $f=gisaidMetadata2Codes($f);
};

#
# vcf2profile
#
$f= vcf2profile($f);

if($opt{out} eq 'profile'){
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
    #   run a bash command taking start stop time
    #
    print "$s";
    print "\nstart: " . qx{date};
    qx{$s};
    print "\nstop: " . qx{date};
}#-----------------------------------

sub vcfListSnippy2Codes{ my ($file) =@_;    
    #-----------------------------------------------
    # Produce a tsv of sample\tvcfcode\tpos 
    # starting from vcf files in snippy format
    # 
    # From
    #   NC_045512	241	snp	C	T	T:2201 C:0	5'UTR	+			intergenic_region n.241C>T 			
    #   NC_045512	344	snp	C	T	T:2269 C:0	mat_peptide	+	79/21290	27/7095	missense_variant c.79C>T p.Leu27Phe	GU280_gp01	ORF1ab	leader protein
    # to
    #   NC_045512      snp:C241T    241
    #   NC_045512      snp:C344T    344
    #-----------------------------------------------
    my $ar=fileList2Array($file);
    my $out= "SAMPLECODE\tVCFCODE\tPOS\n";
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
}#-----------------------------------

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
}#-----------------------------------


sub gisaidMetadata2Codes{ my ($file) =@_; 
    #-------------------------------------
    # Produce a tsv of sample\tvcfcode\tpos 
    # starting from gisaid metadata file format
    # 
    # From
    #   hCoV-19/Tunisia/MHT_2/2020      (N_S202N,NS8_L84S)
    # to
    #   hCoV-19/Tunisia/MHT_2/2020      N_S202N     N_S202
    #   hCoV-19/Tunisia/MHT_2/2020      NS8_L84S    NS8_L84
    # 
    #-------------------------------------
    my $out= "SAMPLECODE\tVCFCODE\tPOS\n";
    open(F, $file);
    while(<F>){
        chomp;
        #print $_;
        ($cmp,$vcfstring)=   _row2_cmp_vcfstring($_);
        # if($_=~/^(\S+)[^\(]+\((\S+)\)/){
        if($cmp){
            #print "$cmp--- $vcfstring\n";
            @vcfs=split(/[,;q]/,$vcfstring);
            my $pos;
            foreach $v (@vcfs) {
                #
                # pos calculation. compatible codes examples
                # GISAID: 
                #   NSP3_I1413L, N_S33del, Spike_ins214EPE  ->  NSP3_1413, N_33, Spike_214
                # INSA algn2pheno script:
                #   N:I1413L,    N:S33del, S:ins214EPE      ->  N:1413,    N:33, S:214
                # NEXTCLADE: 
                #   I1413L,      S33del,   ins214EPE        ->  1413,      33,   214
                #   
                $pos='-1';
                if( $v=~/^\s*([^_:]*[_:]?)(\w*\d+)/ ) {
                    $pos="$1$2";
                }
                #### print "$cmp\t$v\t$pos\n";
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
}#-----------------------------------



sub _row2_cmp_vcfstring{  my ($s) =@_; 
    #
    #   -tsv-sample-pos pos: the position in the tsv file of column containing the sample_name (first position is 0).(default=0)
    #   -tsv-mutationslist-find string=(pos|regexp)*: the way to find the list_of_mutation_codes string in tsv file. if string=pos, -tsv-mutationslist-pos must be set.  if string=regexp, -tsv-mutationslist-regexp must be set. (default=regexp)
    #   -tsv-mutationslist-pos pos*: the position in the tsv file of column containing the list_of_mutation_codes (first position is 0) 
    #   -tsv-mutationslist-regexp string*: the regular expression used to extract the list_of_mutation_codes string. default="`\((.*)\)`"
    #   -tsv-mutation-sep char: the character used as separator between mutations on list_of_mutation_codes string. default=','
    #   -tsv-mutation-pos-regexp string*:  the regular expression used to extract the position of the mutation. default="`^(.*?[_:]?)\w*(\d+)`"
    #   -tsv-mutation-pos-replace string*:
    #
    if( $_=~ /(SAMPLECODE|TYPE|aaSubstitutions|All\smutations)/ ){return ('','');}

    my($cmp,$vcfstring);
    if(                                 $type eq 'gisaid'){
        if($s=~/^(\S+)[^\(]+\((\S+)\)/){
            $cmp=$1;
            $vcfstring=$2;
        }
    }elsif(                             $type eq 'algn2pheno'){
        my @a=split(/\t/,$s);
        $cmp=$a[0];
        $vcfstring=$a[9];

    }elsif(                             $type eq 'nextclade'){
        my @a=split(/\t/,$s);
        $cmp=$a[0];
        $vcfstring=$a[26];
    }elsif(                             $type eq 'tsv' ){
        if( $opt{'tsv-mutationslist-find'} eq 'pos' ){
            my $sep=$opt{'tsv-separator'};
            my $spos=$opt{'tsv-sample-pos'};
            my $mpos=$opt{'tsv-mutationslist-pos'};
            
            my @a=split(/$sep/,$s); 
            $cmp=$a[$spos];
            $vcfstring=$a[$mpos];
        }else{
            my $regexp=$opt{'tsv-mutationslist-regexp'};
            if($s=~/$regexp/){
                $cmp=$1;
                $vcfstring=$2;
            }
        }

    }
    return ($cmp,$vcfstring);
}#-----------------------------------


sub _OLD_gisaidMetadata2Codes{ my ($file) =@_; 
    #-------------------------------------
    # Produce a tsv of sample\tvcfcode\tpos 
    # starting from gisaid metadata file format
    # 
    # From
    #   hCoV-19/Tunisia/MHT_2/2020      (N_S202N,NS8_L84S)
    # to
    #   hCoV-19/Tunisia/MHT_2/2020      N_S202N     N_S202
    #   hCoV-19/Tunisia/MHT_2/2020      NS8_L84S    NS8_L84
    # 
    #-------------------------------------
    my $out= "SAMPLECODE\tVCFCODE\tPOS\n";
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
                # pos calculation. compatible codes examples
                # GISAID: 
                #   NSP3_I1413L, N_S33del, Spike_ins214EPE  ->  NSP3_1413, N_33, Spike_214
                # INSA algn2pheno script:
                #   N:I1413L,    N:S33del, S:ins214EPE      ->  N:1413,    N:33, S:214
                # NEXTCLADE: 
                #   I1413L,      S33del,   ins214EPE        ->  1413,      33,   214
                #   
                $pos='-1';
                if( $v=~/^(.*?[_:]?)\w*(\d+)/ ) {
                    $pos="$1$2";
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
}#-----------------------------------

#############################################
#  vcf -> profile
#############################################


sub vcf2profile { my ($f)=@_;
    #-------------------------------------
    #   Hamming Distance matrix from Samples vcf codes
    #   usage: 
    #   vcf2mst.pl samples_vcfcodes.csv > hamming_distance_matrix.tsv
    #
    #   print a matrix compatible with grapetree input
    #-------------------------------------
    my $out_file="/tmp/hamming_distance_matrix.tsv";

    my $c={};
    my $codes={};
    my @a; 

    open(F, $f);
    while(<F>){
        chomp;
        # check headers
        if( $_=~ /(SAMPLECODE|TYPE)/ ){next;}
        #
        # extract code, mutations and pos
        #
        my ($SAMPLECODE,$VCFCODE,$POS) =   split(/\t/,$_);
        ###########àprint "($SAMPLECODE,$VCFCODE,$POS)\n";# =   split(/\t/,$_);
        
        if( $VCFCODE eq ''){
            ($SAMPLECODE,$VCFCODE,$POS) =   split(/,/,$_);
        }        
        if( ! exists($c->{$SAMPLECODE})){
            $c->{$SAMPLECODE}={};
        }
        #
        # filter mutation  based on position 
        #
        if( _avoid_mutation_by_pos($POS) ){
            next;
        }
        # 
        # "allele" code based on variant position
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
    }
    close(F);
    #---------------------------
    # out
    #---------------------------
    my $out="";
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
    #
    # print
    #
    open(F, ">$out_file" );
    print F $out;
    close(F);    
    return $out_file;
}#-----------------------------------

sub _avoid_mutation_by_pos{ my ($pos)=@_;
    # TODO: manage -minmax and -minmaxExclude
    return ''; #false
}#-----------------------------------



sub init{
    #--------------------------------------
    #   INIT
    #--------------------------------------
    if( ! $f ){
        print $presentation;
        exit;          
    }    
    if( $type =~ /^(vcf|gisaid|algn2pheno|nextclade|tsv|code)$/ ){
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

    #
    # options
    # 
    init_set_option_defaults();
    my $initial_options=  join(' ',@ARGV);
    $initial_options   =~ s/^.*?-/-/;
    ## print "initial_options=$initial_options\n";
    $options=$initial_options;
    while( $options =~ /^\s*-(\S+)\s+(\S+)(.*)$/ ){
        print "---options=$options\n";
        $opt{$1}=$2;
        $options=$3;
    }
    if( $initial_options ne $options){
        print "find options:\n";
        foreach my $k (keys(%opt)){

            if(   $k eq 'out' ){
                if($opt{$k} !~ /profile/){
                    print qq{ERR: value "$opt{$k}" not known for option "$k"
                    example: -out profile\n};
                    exit;
                }
            }
            elsif($k =~ /^(minmax|minmax-exclude)$/ ){
                if($opt{$k} !~ /^([^:,]+:[^:,]+,?)+$/){
                    print qq{ERR: value "$opt{$k}" not correct for option "$k"
                    example: -minmax 0-1000
                    example: -minmax 0-1000,1200-12111
                    \n};
                    exit;
                }
            }
            elsif($k =~ /^(tsv-separator|tsv-sample-pos|tsv-mutationslist-find|tsv-mutationslist-pos|tsv-mutationslist-regexp|tsv-mutation-sep|tsv-mutation-pos-regexp|tsv-mutation-pos-replace)$/){
                #ok
            }else{
                print qq{ERR: option "$k" not recognized \n};
                exit;
            }
        }
    }
   
    if($options){
        print "options not recognized\n$options\n";
        exit;
    }
}


sub init_set_option_defaults {
    $opt{'tsv-separator'}="\t";
    $opt{'tsv-sample-pos'}=0;
    $opt{'tsv-mutationslist-find'}='regexp';
    $opt{'tsv-mutationslist-pos'}=0;
    $opt{'tsv-mutationslist-regexp'}='\((.*)\)';
    $opt{'tsv-mutation-sep'}=',';
    $opt{'tsv-mutation-pos-regexp'} ='^(.*?[_:]?)\w*(\d+)';
    $opt{'tsv-mutation-pos-replace'}='$1$2';
};

<<_________COMMENT_____________;


* *-minmax value*: Take mutations with position in the "`value`" interval. Format `value` is `min1:max1,min2:max2`. Example `-minmax 0-100,200-1500,5000-5500`
* *-minmaxExclude value*: Exlude mutations with position in the "`value`" interval. Format is the same of `-minmax` value

* *-tsv-XXX*: different options for manipulating a tsv file containing at least 2 columns sample_name and list_of_mutation_codes. This options are considered only in case of *type_of_input=tsv*
  * *-tsv-separator char: the character used as separator on tsv/csv file.(default='\t')
  * *-tsv-sample-pos pos: the position in the tsv file of column containing the sample_name (first position is 0).(default=0)
  * *-tsv-mutationslist-find string=(pos|regexp)*: the way to find the list_of_mutation_codes string in tsv file. if string=pos, -tsv-mutationslist-pos must be set.  if string=regexp, -tsv-mutationslist-regexp must be set. (default=regexp)
  * *-tsv-mutationslist-pos pos*: the position in the tsv file of column containing the list_of_mutation_codes (first position is 0) 
  * *-tsv-mutationslist-regexp string*: the regular expression used to extract the list_of_mutation_codes string. default="`\((.*)\)`"
  * *-tsv-mutation-sep char: the character used as separator between mutations on list_of_mutation_codes string. default=','
  * *-tsv-mutation-pos-regexp string*:  the regular expression used to extract the position of the mutation. default="`^(.*?[_:]?)\w*(\d+)`"
  * *-tsv-mutation-pos-replace string*: the regular expression used to extract the position of the mutation. default="`$1$2`"



            if($k =~ /^(minNT|max-nt|min-aa|max-aa)xNT
           
                print " -$k=$opt{$k}\t{minNT};{
                    if($POS ><$k=$o ){return;}
next }
            }elsprint " -$k=$opt{$k}\t{maxNT};{
                    if($POS ><$k=$o ){retur>;}
next          maxNT
            }else{
                print " -$k is not recgnized\n";
            }
        }
    }  
algn 9
next 26

_________COMMENT_____________


