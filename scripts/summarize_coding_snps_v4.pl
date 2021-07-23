#!/usr/bin/perl
use strict;

my $usage = "$0 -i <vcf file> -f <alt_allele_freq [0.8]> -p T/F [print only high quality alleles] -m T/F [print missense alleles only] -n <gene_name_file>\n\n";
my %arg = @ARGV;
die $usage unless ($arg{-i} and $arg{-p} and $arg{-m});
my %names;
my $ALLELE_FREQUENCY = $arg{-f} || 0.8; ## minimum allele frequency to call variants in each sample
if ($arg{-n}){
	open (NAME,"<$arg{-n}");
	while(<NAME>){
		chomp;
		my @x=split /\t/;
		$names{$x[0]} = $x[1];	
	}
	close NAME;
}

open (FHI, "<$arg{-i}") || die "ERROR, I cannot open $arg{-i}: $!\n\n";
#print "CHROMOSOME\tPOSITION\tREFERENCE\tVARIANT\tCOVERAGE\tALLELE_FREQ\tQUAL_FILTER\tMUTATION_TYPE\tMUTATION\tGENE_ID\tPRODUCT_NAME\n";
my (@freq, @depth);
while(<FHI>){
	chomp;
	if(m/^#CHROM/){
		my @fields = split /\t/;
		my $l = scalar(@fields)-1;
		print "CHROMOSOME\tPOSITION\tREFERENCE\tVARIANT\t";
		foreach my $sample (@fields[9..$l]){
			print "COVERAGE_$sample\tALLELE_FREQ_$sample\tALLELE_$sample\t";
		}
		print "QUAL_FILTER\tMUTATION_TYPE\tMUTATION\tGENE_ID\tPRODUCT_NAME\n";
	}
	next if m/^#/;
	my @x = split /\t/;
	my $depth_freq;
	
	my ($chrom, $pos, $ref, $alt) = ($x[0], $x[1],$x[3], $x[4]);

	## Number of columns?
	my $l = scalar(@x)-1;
	foreach my $sample (@x[9..$l]){
	 
                my ($freq, $depth, $genotype) = &get_freq($sample);
                my $allele;
                if($PLOIDY == 1){
                        #$allele = $freq >= $ALLELE_FREQUENCY ? $alt : $ref;
                        $allele = $genotype  == 1? $alt : $ref;
                }
                if($PLOIDY == 2){
                        $allele = $freq <= $MIN ? "0_0" : $freq > $MAX ? "1_1" : "0_1";
                }
                $depth_freq .= $depth."\t".$freq."\t$allele\t";
	}
	
	my $pass = $x[6];
	my ($pub_locus, $mutation, $type) = &get_mutation($x[7]);
	#my ($chrom, $pos, $ref, $alt) = ($x[0], $x[1],$x[3], $x[4]);
#	print "pass=$pass freq=$freq type=$type  pub_locus=$pub_locus, mutation=$mutation\n";
	next if ( ($arg{-p} eq 'T') and ($pass ne 'PASS') );
	next if ( ($arg{-m} eq 'T') and !( ($type =~m/missense_variant/) or ($type =~ m/splice/ ) or ($type =~ m/stop/ ) ) );
	# next if ( $freq < $arg{-f} );
	
	## print output
	print "$chrom\t$pos\t$ref\t$alt\t$depth_freq"."$x[6]\t$type\t$mutation \t$pub_locus\t$names{$pub_locus}\n";
}
close FHI;
exit(0);

##############################
sub get_freq {
        my $data = shift;
        my @x = split(':',$data);
        my $AD = $x[1];
        my ($AD_R, $AD_A) = split(',' , $AD);
        my $DP = $x[2];
        my $GT = $x[0];
        #my ($r,$a)= ($x[2], $x[3]); # ($1, $2) if $x[1]=~m/^(\d+),(\d+)/;
        my $total = $DP ;
        #print "$data\n" if $total == 0;;
        my $freq = $total ? $AD_A / $DP : 0;
        #print "a = $a , r = $r , $data\n";
        $freq = $1 if $freq =~ m/^(\d+\.?\d?\d?)/;
        return ($freq, $DP, $GT);
}

sub get_mutation {
	my $data = shift;
	my @stats = split(';', $data);
	my @annot = split(m/\|/, $data);
	#print "$stats[-1]\n";
	my ($pub_locus, $nuc_mutation, $mutation, $type) = ($annot[3], $annot[9],  $annot[10], $annot[1]);
	$pub_locus =~ s/\s+//g;
	my $my_mutation = $mutation || $nuc_mutation;
	return ($pub_locus, $my_mutation, $type) ;
}

# INDEL;IDV=24;IMF=0.648649;DP=121;VDB=4.94589e-05;SGB=-2.73308;MQSB=0.120735;MQ0F=0;DPR=0,106;AF1=1;AC1=8;DP4=0,0,13,93;MQ=17;FQ=-90.9422;ANN=TGGG|frameshift_variant|HIGH|TGGT1_220650|TGGT1_220650|transcript|rna_TGGT1_220650-1|Coding|5/5|c.1033dupG|p.Glu345fs|1034/1110|1034/1110|345/369||WARNING_TRANSCRIPT_NO_START_CODON;LOF=(TGGT1_220650|TGGT1_220650|1|1.00)

