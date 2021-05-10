#-------------------------------------
#   vcflist2codes.pl 
#  
#   Hamming Distance based Minimum Spanning Tree from Samples vcf
#   usage 1: 
#   vcflist2codes.pl list_of_vcfiles  > samples_vcfcodes.csv
#

#   list_of_vcfiles = list of files. 
#   Each file in the list must have the vcf format like this (e.g., from snippy output):

#   #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  DS11243446_vdsnippy_NC045512
#   NC_045  241     .       C       T       42316.5 .       AB=0;AO=1409;DP=1411;QA=47794;QR=36;RO=1;TYPE=snp;ANN=T|intergenic_region|MODIFIER|CHR_START-GU280_gp01|CHR_START-GENE_GU280_gp01|intergenic_region|CHR_START-GENE_GU280_gp01|||n.241C>T||||||  GT:DP:RO:QR:AO:QA:GL    1/1:1411:1:36:1409:47794:-47794
#-------------------------------------


my ($f)=@ARGV;

while(<>){
    chomp;
    $f=$_;
    $cmp=$f;

    print "SAMPLECODE,VCFCODE\n";

    open(F, $f);
    while(<F>){
        if( $_ =~ /^#/){next;}

        my $r=$_;
        my @ar=split(/\t/,$r);
        my $type='';
        if( $r =~ /;TYPE=(\w+);/){
            $type=$1;
        }


        if( $type eq 'snp'){
            print "$cmp,$type:@ar[3]@ar[1]@ar[4]\n";
        }elsif( $type eq 'mnp'){
            print "$cmp,$type:@ar[3]@ar[1]@ar[4]\n";
        }elsif( $type =~ 'del' ){
            $pos=@ar[1] + 1;
            $pos2=$pos + length(@ar[3]) -2;
            print "$cmp,$type:$pos-$pos2\n";
        }else{
            print "$cmp,$type:@ar[1]-@ar[3]-@ar[4]\n";
        }
    }
    close(F);
}

