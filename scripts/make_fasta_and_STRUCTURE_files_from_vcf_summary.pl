#!/usr/local/bin/perl
use strict;
use Data::Dumper;

my $usage = "$0 -p <output file prefix> -R <reference_name/F [F] include reference column as reference> -i <vcf summary file> -c <minimum coverage [5]> -min <minimum alt allele freq range [1]> -max <maximum alt allele freq range [0]> -M <min average coverage to include sample [0]> -d <discard intervals chr:end5:end3,chr:end5:end3,...>\n\n";
my %arg = @ARGV;
die $usage unless $arg{-i} and $arg{-p};
my $COV = $arg{-c} || 5;
my $MIN = $arg{-min} || 1;
my $MAX = $arg{-max} || 0;
my $REF = $arg{-R} || 0; #? $arg{-R} : 0;
die "ERROR, minimum and maximum alt allele frequency value should be between 0 and 1. For example 0.8 for 80%\n\n" if $MIN > 1 or $MAX > 1;
my $int = $arg{-d} || 0;

my %int;
my $MIN_COV = $arg{-M} || 0;
my $fasta = $arg{-p}.".fasta";
my $structure = $arg{-p}.".structure.txt";

if(-e $fasta){ die "ERROR: file $fasta already exists\n\n" }elsif(-e $structure){die "ERROR: file $structure already exists\n\n"}

open (FASTA,  ">$fasta");
open (STR, ">$structure");

my @int = split("," , $arg{-d});
foreach my $i (@int){
    my @feat = split(':' , $i);
    push @{ $int{$feat[0]} }, $feat[1]."_".$feat[2]; 
    #print ">> @feat\n";
}

open (SNP, "<$arg{-i}") || die "ERROR, I cannot open $arg{-i}: $!\n\n";
my $header = <SNP>; chomp($header);
#my @x = split(/\t/,$header);
my @samples;
while($header =~ m/COVERAGE_(\S+)/g){
    push @samples, $1;
}

my $num_snps;
my %cov;
my @snps;
while(<SNP>){
    push @snps, $_;
    $num_snps++;
    my $c = 0;
    my @x = split /\t/;
    $cov{$REF} += 1;
    foreach my $sample (@samples){
        $cov{$sample} += $x[4+($c*3)];
        $c++;
    }
}
close SNP;

## Calculate average cov
foreach my $sample (@samples){
    $cov{$sample} = $cov{$sample} / $num_snps;
}
$cov{$REF} = $cov{$REF} / $num_snps if $REF;

my %seq; ## holds sequences
my @ids; ## hods SNP ids in same order as SNPs stored in %seq
foreach my $snp (@snps){
    my @x = split(/\t/, $snp);
    my $chr = shift @x;
    my $pos = shift @x;
    my $ref = shift @x;
    my $alt = shift @x;
    my $SNP_ID = $chr."_".$pos;
    #print ">$chr - $pos\n";
    ## skip InDels
    next if (length($ref) > 1 or length($alt) > 1 or $ref eq '-' or $alt eq '-');
    
    ## skip SNPs within discarded intervals from $arg{-d}
    my $flag_excl = 0;
    #print ">>> $int{$chr}\n";
    foreach my $pair ( @{ $int{$chr} } ){
        my ($end5, $end3) = split("_", $pair);
        $flag_excl = 1 if ($pos >= $end5 and $pos <= $end3);
        #print "$pos >= $end5 and $pos <= $end3 \n";
    }
    next if $flag_excl == 1;

    my $flag = 0 ; ## set $flag = 1 if some of the coverage of allele freq conditions are not fulfilled
    my %seqTmp; # holds temporary sequence data for a given position
    $seqTmp{$REF} = $ref if $REF;
    foreach my $sample (@samples){
        my $cov = shift @x;
        my $freq = shift @x;
        my $allele = shift @x;
        next if $cov{$sample} < $MIN_COV;
        $flag = 1 if $cov < $COV;
        $flag = 1 if ( $freq > $MIN and $freq < $MAX );
        $seqTmp{$sample} = $allele;
    }
    next if $flag == 1; ## conditions not fulfilled

    #	print 	Dumper(%seqTmp);;
    
    ## add temp seq if everything looks fine
    foreach my $sample (@samples){
        $seq{$sample} .= $seqTmp{$sample};
    }
    $seq{$REF} .= $seqTmp{$REF} if $REF;
    push @ids, $SNP_ID;
}
close SNP;


## Print out fasta file
@samples = ($REF, @samples);
foreach my $sample (@samples){
    next if $cov{$sample} < $MIN_COV;
    print FASTA ">$sample\n$seq{$sample}\n";
}

## Print structure file
print STR "\t".join("\t", @ids)."\n";
foreach (@ids){print STR "\t0.01"};print STR "\n";
foreach my $sample (@samples){
    next if $cov{$sample} < $MIN_COV;
    print STR "$sample";
    foreach my $snp (split(m//, $seq{$sample}) ){
        print STR "\t$snp";
    }
    print STR "\n";
}

close FASTA;
close STR;

