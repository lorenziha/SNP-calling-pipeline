#!/usr/local/bin/perl
use strict;

######################################################################################################
#
#  The pipeline assumes 32 CPUs available
#
######################################################################################################

# SET THE FOLLOWING PATHS
my $trimmomatic = "/usr/local/devel/ANNOTATION/hlorenzi/bin/trimmomatic.jar";
my $bowtie = "/usr/local/devel/ANNOTATION/hlorenzi/bin/bowtie2";
my $bwa = "/usr/local/bin/bwa";
my $sortsam = "/usr/local/devel/ANNOTATION/hlorenzi/picard-master/dist/picard.jar SortSam";
my $markduplicates = "/usr/local/devel/ANNOTATION/hlorenzi/picard-master/dist/picard.jar MarkDuplicates";
my $samtools = "/usr/local/bin/samtools";
my $GATK = "/usr/local/devel/ANNOTATION/hlorenzi/bin/GenomeAnalysisTK.jar";
my $java2 = "/usr/local/devel/ANNOTATION/hlorenzi/jre1.7.0_75//bin/java";
my $java = "/usr/local/packages/jdk-8u66/bin/java"; #"/usr/local/bin/java";
my $bcftools = "/usr/local/bin/bcftools";
my $filter = "-sLowQual -g3 -G10 -e'%QUAL<10 || (%MAX(DV)<=3 && FMT/GT=\"./1\")|| (%MAX(DV)/%MAX(DP)<=0.3 && FMT/GT=\"./1\") || %MAX(FMT/DP)<5 || (RPB>=0 && RPB<0.1 && %QUAL<15) || (MQB>=0 && MQB<0.05) || (BQB>=0 && BQB<0.05) || (MQSB>=0 && MQSB<0.0)'";
my $snpEff = "/usr/local/devel/ANNOTATION/hlorenzi/snpEff/snpEff.jar";
my $summarize_coding_snps = "/usr/local/devel/ANNOTATION/hlorenzi/PERL/SCRIPTS/summarize_coding_snps_v4.pl";
my $adapters = "/usr/local/devel/ANNOTATION/hlorenzi/Trimmomatic-0.33/adapters/adapters_pipeline";
my %log; 
my $bouwtiw2build = "/usr/local/devel/ANNOTATION/hlorenzi/bowtie2/bowtie2-build";
my $picard = "/usr/local/devel/ANNOTATION/hlorenzi/picard-master/dist/picard.jar";
##################################################################################################


my $usage = "$0 -r <reference_genome> [-nosnp T/F default F][-o <oligos_file>] [-R <parent/reference fastq prefix if any>] -S [skip mappings, just do SNP calls, defaulf = F]  -D T/F [delete temp files default F] -f <file with other fastq file prefixes> [-c chromosome_ID] -A <reference annotation [RH88/GT1/ME49/Pf3D7/ATH/BESB1 other path to file]> \n\n";

# Check options
my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
my %arg = @ARGV;
die $usage unless ($arg{-r} && $arg{-f} && $arg{-A});
my $NOSNP = $arg{-nosnp} eq 'T'? 1 : 0 ;
my $DELETE = $arg{-D} eq 'T'? 1:0;
my $SKIP_MAPPING = $arg{-S} eq 'T'? 1:0;

my $CWD = `pwd`; chomp $CWD;

my ($chr, $L);
$chr = $arg{-c}.'.' if $arg{-c} ;
$L = "-L $arg{-c} " if $arg{-c};


## CHECK LOG FOR PIPELINE STATUS (RESTART OR RESUME)
if (-e "snp_call.log"){
	open (LOG, "<snp_call.log");
	while(<LOG>){
		if(m/^STEP:(\S+)/){
			$log{$1} = 1;
		}
	}

}

## Use custom SNP filtering parameters
$filter = $arg{-F} if $arg{-F};


## Tabulated file: Gene_ID <TAB> Product_name
## Gene_ID should match Gene_ID in gff3 file used for snpEff

## Available genome annotation files (accesion number <TAB> Description/Gene name)
my %annotation = ( "RH88" => "/usr/local/devel/ANNOTATION/hlorenzi/snpEff/data/Tgondii_RH88/rh88_com_name.txt",	
		   "GT1" => "/usr/local/projdata/700030/projects/GCID_PARASITES/KATIE_CLONES/tggt1_names.txt",
		   "ME49" => "/usr/local/devel/ANNOTATION/hlorenzi/snpEff/data/ToxoDB-13.0_TgondiiME49/tga4_com_name",
		   "Pf3D7" => "/usr/local/devel/ANNOTATION/hlorenzi/snpEff_v4_1b/data/PlasmoDB_30_Pfalciparum3D7/Pf3D7_names.txt",
		   "ATH" => "/usr/local/devel/ANNOTATION/hlorenzi/snpEff_v4_1b/data/athalianaTair10/athalianaTair10_names.txt",
		   "BESB1" => "/usr/local/devel/ANNOTATION/hlorenzi/snpEff_v4_1b/data/besb1.11012017.assembly/besb1_names.txt",
	           "BS" => "/usr/local/devel/ANNOTATION/hlorenzi/snpEff_v4.3/data/BS168/BS168_com_name.txt"
		);
my $annotation_file = $annotation{ $arg{-A} }? $annotation{ $arg{-A} } : $arg{-A};
die "ERROR, I cannot find file $arg{-A}\n" unless (-e $annotation_file);

## Load prefixes
my @prefix = ($arg{-R}) || ();
if ($arg{-f}){
	open (PREFIX,"<$arg{-f}") || die "ERROR, I cannot open $arg{-f}: $!\n\n";
	while(<PREFIX>){
		chomp;
		next if m/^[\s\t]*$/;
		s/\s+//g;
		push  @prefix, $_;
	}
}



my $gz;
my $bwa_ref_prefix_orig = $arg{-r};
my $bwa_ref_prefix = $bwa_ref_prefix_orig;
$bwa_ref_prefix =~ s/\.(fasta|fa|fsa)$//;

## Chek bwa index files for. reference are there, otherwise create them:
unless (-e $bwa_ref_prefix_orig.'.sa'){
        warn "I cannot find index files for reference genome $arg{-r}.\n"
        ."Creating Index files with base name: $bwa_ref_prefix\n" unless $SKIP_MAPPING;
        my $CMD = "$bwa index  $arg{-r}";
        &print_stderr ("## Building bwa index files with prefix $bwa_ref_prefix",$CMD) unless $SKIP_MAPPING;
}

## Chek dictionary files are there, otherwise create them:
unless (-e $bwa_ref_prefix.'.dict'){
	my $CMD = "$java -jar $picard CreateSequenceDictionary REFERENCE=$arg{-r} O=$bwa_ref_prefix.dict";
	&print_stderr ("## Building dict file $bwa_ref_prefix.dict",$CMD) ;
}

## Chek fai files are there, otherwise create them:
unless (-e "$arg{-r}.fai"){
	my $CMD = "samtools faidx $arg{-r}";
	 &print_stderr ("## Building fai file $arg{-r}.fai",$CMD);
	my $CMD = "ln -s $arg{-r}.fai $bwa_ref_prefix.fai";
	 &print_stderr ("## Linking $arg{-r}.fai to $bwa_ref_prefix.fai",$CMD) ;
}


## Chek dict index files are there, otherwise create them::
unless (-e $bwa_ref_prefix.'.dict'){
       die "ERROR, I cannot find dictionary dict files for reference genome $arg{-r}.\nIndex files name should be: $bwa_ref_prefix.dict\n";
};
## Chek fai index files are there:
unless (-e "$arg{-r}.fai"){
       die "ERROR, I cannot find fai index files for reference genome $arg{-r}\nIndex files name should be: $bwa_ref_prefix.fai\n";
 };

## Trim fastq files
foreach my $prefix (@prefix){
	next if $SKIP_MAPPING;
	my $flag = "trimmomatic.".$prefix;
	$gz = (-e $prefix.".R2.fastq.gz") ? "fastq.gz" : "fastq";

	## next if (-e $prefix.".paired.R2.$gz"); ## skip file if trimmed file already exists DEPRECATED
	my $CMD = "$java -jar $trimmomatic PE -threads 32 -phred33 $prefix.R1.$gz $prefix.R2.$gz $prefix.paired.R1.$gz $prefix.unpaired.R1.$gz $prefix.paired.R2.$gz $prefix.unpaired.R2.$gz ILLUMINACLIP:$adapters:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50 LEADING:12 TRAILING:3 "; 

	#print STDERR "## Trimming reads from $prefix.$gz\n$CMD\n\n";
	&print_stderr ("## Trimming reads from $prefix.$gz file",$CMD) unless $log{$flag};
	&print_log($flag) unless $log{$flag};
	if($DELETE){
		`rm -f "$prefix.R1.$gz"`;
		`rm -f "$prefix.R2.$gz"`;
	}

}

 
## Map reads to reference and remove unmapped reads from bam file
foreach my $prefix (@prefix){
	next if $SKIP_MAPPING;
	my $flag = "bwa_mem.".$prefix;
	$gz = (-e $prefix.".paired.R2.fastq.gz") ? "fastq.gz" : "fastq";
	my $CMD = "$bwa mem -R ".'"@RG\tID:'.$prefix.'\tPL:Illumina\tSM:'.$prefix.'\tLB:'.$prefix.'" '." -t 32 $arg{-r} $prefix.paired.R1.$gz $prefix.paired.R2.$gz | $samtools view -bSu -F 4 - > $prefix.bam";
	&print_stderr ("## Mapping reads from $prefix to reference",$CMD) unless $log{$flag};
	&print_log($flag) unless $log{$flag};
	## I am not deleting the trimmed reads here

	# remove unmapped reads from bam
	my $CMD = "samtools view -@ 32 -h -b -F 4 $prefix.bam > $prefix.mapped.bam";
	my $flag = "remove_unmapped_reads.$prefix";
	&print_stderr ("## Removing unmapped reads from $prefix.bam to reference",$CMD) unless $log{$flag};
	&print_log($flag) unless $log{$flag};

	# generate files with unmapped reads
	my $CMD = "samtools view -@ 32 -h -b -f 4 $prefix.bam > $prefix.unmapped.bam";
	my $flag = "create_unmapped_read_file.$prefix";
	&print_stderr ("## Removing unmapped reads from $prefix.bam to reference",$CMD) unless $log{$flag};
	&print_log($flag) unless $log{$flag};
}

## Mark duplicates from bam
my @files_to_map;
foreach my $prefix (@prefix){
	my $flag = "sortsam.".$prefix;
	push @files_to_map, "$prefix.raw.snps.indels.g.vcf";
	next if $SKIP_MAPPING;

	## Sort bam file
	my $CMD = "samtools sort -O bam -@ 32 -T tmp.bam -o $prefix.sorted.bam $prefix.mapped.bam"; 
	&print_stderr ("## Sorting $prefix.bam file",$CMD) unless $log{$flag};
	&print_log($flag) unless $log{$flag};
	if($DELETE){
		`rm -f "$prefix.bam"`;
		`rm -f "$prefix.mapped.bam"`;
	}
	
	## Mark duplicated reads
	my $flag = "dedup.".$prefix;
	my $CMD = "$java -jar $markduplicates INPUT=$prefix.sorted.bam OUTPUT=$prefix.sorted.dedup.bam METRICS_FILE=$prefix.sorted.metrics  REMOVE_DUPLICATES=true";
	&print_stderr ("## Markind duplicates in $prefix.sorted.bam file",$CMD) unless $log{$flag};
	&print_log($flag) unless $log{$flag};
	if($DELETE){
		`rm -f "$prefix.sorted.bam"`;
		`rm -f "$prefix.sorted.bai"`;
		`rm -f "$prefix.sorted.metrics"`;
	}	
	
	## indexing
	my $flag = "indexing.".$prefix; 
	my $CMD = "$samtools index $prefix.sorted.dedup.bam";
	&print_stderr ("## indexing $prefix.sorted.dedup.bam file",$CMD) unless $log{$flag};
	&print_log($flag) unless $log{$flag};
	## I am not deleting $prefix.sorted.dedup.bam file here

	if($DELETE){
		`rm -f "$prefix.sorted.dedup.bam"`;
		`rm -f "$prefix.sorted.dedup.bai"`;
	}
}


## STOP pipeline here if no SNP call required
exit(0) if $NOSNP;

## Look for SNPs with haplotype calling
foreach my $prefix (@prefix){
	# haplotype calling
	my $flag = "haplotypecaller.".$prefix;
	my $CMD = "$java -jar $GATK -T HaplotypeCaller -R $arg{-r} -I $CWD/$prefix.sorted.dedup.bam --emitRefConfidence GVCF -ploidy 1 --dontUseSoftClippedBases --pcr_indel_model NONE  -o $CWD/$prefix.raw.snps.indels.g.vcf -XL $CWD/Tgondii_RH88.exclude.intervals -XL $CWD/unitig.intervals --num_cpu_threads_per_data_thread 30";
	&print_stderr ("## Running HaplotypeCaller on $prefix.sorted.dedup.bam file",$CMD) unless $log{$flag};
	&print_log($flag) unless $log{$flag};

}

## SNP calls with GenotypeGVCFs
my $p = join("_",@prefix);

## check if file name will be too long
$p = length($p) > 50 ? "ALL_SAMPLES" : $p;
$p = $chr.$p;
my $flag = "GenotypeGVCFs.$chr".$p;
my @gfiles_to_merge;
foreach my $file (@files_to_map){
	push @gfiles_to_merge, ("--variant", $file);
}
my $CMD = "$java -Xmx200G  -jar $GATK -T GenotypeGVCFs -R $arg{-r} -nt 30 -o $p.raw.snps.indel.vcf $L -XL Tgondii_RH88.exclude.intervals ".join(" ",@gfiles_to_merge);
&print_stderr ("## Running GenotypeGVCFs on @gfiles_to_merge",$CMD) unless $log{$flag};
&print_log($flag) unless $log{$flag};

## SNP filtering
$flag = "SelectVariants_SNP".$p;
my $CMD = "$java -Xmx200G  -jar $GATK -T SelectVariants -R $arg{-r} -V $p.raw.snps.indel.vcf -selectType SNP -o $p.raw_snps.vcf";
&print_stderr ("## Running SNP extraction on $p.raw.snps.indel.vcf",$CMD) unless $log{$flag};
&print_log($flag) unless $log{$flag};

## NOTE: I eliminated ReadPosRankSum and MQRankSum from the filter because those tests require at least one sample to be heterozygote and that is not possible in most parasite genomes like Apicomplexans

$flag = "VariantFiltration_SNP".$p;
my $CMD = "$java -Xmx200G -jar $GATK -T VariantFiltration -R $arg{-r} -V $p.raw_snps.vcf --filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 \" --filterName \"my_snp_filter\" -o $p.filtered.snps.vcf";
#my $CMD = "$java -jar $GATK -T VariantFiltration -R $arg{-r} -V $p.raw_snps.vcf --filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum <     -12.5 || ReadPosRankSum < -8.0\" --filterName \"my_snp_filter\" -o $p.filtered.snps.vcf";
&print_stderr ("## Running SNP filtering  on $p.raw_snps.vcf",$CMD) unless $log{$flag};
&print_log($flag) unless $log{$flag};

## InDels filtering
$flag = "SelectVariants_InDel".$p;
my $CMD = "$java -Xmx200G  -jar $GATK -T SelectVariants -R $arg{-r} -V $p.raw_snps.vcf -selectType INDEL -o $p.raw_indels.vcf";
&print_stderr ("## Running InDel extraction  on $p.raw.snps.indel.vcf",$CMD) unless $log{$flag};
&print_log($flag) unless $log{$flag};

$flag = "VariantFiltration_InDel".$p;
my $CMD = "$java -Xmx200G -jar $GATK -T VariantFiltration -R $arg{-r} -V $p.raw_indels.vcf --filterExpression \"QD < 2.0 || FS > 200.0 \" --filterName \"my_indel_filter\" -o $p.filtered.indels.vcf";
#my $CMD = "$java -jar $GATK -T VariantFiltration -R $arg{-r} -V $p.raw_indels.vcf --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.    0\" --filterName \"my_indel_filter\" -o $p.filtered.indels.vcf";
&print_stderr ("## Running InDel filtering  on $p.raw_snps.vcf",$CMD) unless $log{$flag};
&print_log($flag) unless $log{$flag};

## Generate snpEff-based annotation for SNPs and INDELS separated
my $flag = "snpEff_SNP".$p;
my $ref = $arg{-A} eq "RH88" ? "Tgondii_RH88" : $arg{-A} eq "GT1" ? "gt1_genome_w_RH_api" : $arg{-A} eq "ME49" ? "ToxoDB-13.0_TgondiiME49" : $arg{-A} eq "Pf3D7" ? "PlasmoDB_30_Pfalciparum3D7" :  $arg{-A} eq "ATH" ? "athalianaTair10" : $arg{-A} eq "BESB1" ? "besb1.11012017.assembly" : "OTHER";
if ($ref eq "OTHER" ){
	print STDERR "No available library for $arg{-A}. SnpEff will not run\n\n";
	exit(0);
}
my $CMD = "$java -Xmx4g -jar $snpEff -stats $p.snpEff_summary.snps.html -ud 1000 $ref $p.filtered.snps.vcf  > $p.filtered.snps.annotation.vcf";
&print_stderr ("## Running snpEff on $p.filtered.snps.vcf",$CMD) unless $log{$flag};
&print_log($flag) unless $log{$flag};

my $flag = "snpEff_INDELS".$p;
my $CMD = "$java -Xmx4g -jar $snpEff -stats $p.snpEff_summary.indels.html -ud 1000 $ref $p.filtered.indels.vcf  > $p.filtered.indels.annotation.vcf";
&print_stderr ("## Running snpEff on $p.filtered.indels.vcf",$CMD) unless $log{$flag};
&print_log($flag) unless $log{$flag};

my $flag = "summarize_coding_snps".$p;
my $CMD = "perl $summarize_coding_snps -i $p.filtered.snps.annotation.vcf -p F -m F -n $annotation_file > $p.filtered.snps.annotation.summary.txt";
&print_stderr ("## Running summarize_coding_snps.pl on $p.filtered.snps.annotation.vcf",$CMD) unless $log{$flag};
&print_log($flag) unless $log{$flag};

my $flag = "summarize_coding_indels".$p;
my $CMD = "perl $summarize_coding_snps -i $p.filtered.indels.annotation.vcf -p F -m F -n $annotation_file > $p.filtered.indels.annotation.summary.txt";
&print_stderr ("## Running summarize_coding_snps.pl on $p.filtered.indels.annotation.vcf",$CMD) unless $log{$flag};
&print_log($flag);

exit(0);

####################################
# Functions
####################################

sub print_stderr {
	my ($message,$cmd) = @_;
	print STDERR "\n\n","="x100,"\n$message\n$cmd\n","="x100,"\n\n";
	my $err = system("$cmd");
	#`$cmd`;
	if ($err){
		die "ERROR, the commands below failed with error $err\n$cmd\n\n";
	}
}

sub print_log {
	my $flag = shift;
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
	open (LOG, ">>snp_call.log") || die "ERROR, I cannot open snp_call.log file: $!\n\n";
	print LOG "STEP:".$flag." $hour:$min $mday $months[$mon] $days[$wday]\n";
	close LOG;
}
