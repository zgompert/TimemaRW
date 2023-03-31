#!/usr/bin/perl
#
# conver sam to bam, then sort and index 
#


use Parallel::ForkManager;
my $max = 24;
my $pm = Parallel::ForkManager->new($max);

## green striped 
#my $genome = "/uufs/chpc.utah.edu/common/home/u6000989/data/timema/hic_genomes/t_crist/timema_cristinae_12Jun2019_lu3Hs.fasta";
## green unstriped
my $genome = "/uufs/chpc.utah.edu/common/home/gompert-group1/data/timema/hic_genomes/t_crist_gus/mod_hic_output.fasta";

FILES:
foreach $fq (@ARGV){
	$pm->start and next FILES; ## fork
        if ($fq =~ m/(2013[A-Z0-9_]+)/){
        	$ind = $1;
    	}
    	else {
       		die "Failed to match $file\n";
    	}
	system "gunzip $fq\n";
	$fq =~ s/\.gz// or die "failed sub\n";
	system "bwa aln -n 4 -l 20 -k 2 -t 1 -q 10 -f aln"."$ind".".sai $genome $fq\n";
	system "bwa samse -n 1 -r \'\@RG\\tID:tcr-"."$ind\\tPL:ILLUMINA\\tLB:tcr-"."$ind\\tSM:tcr-"."$ind"."\' -f aln"."$ind".".sam $genome aln"."$ind".".sai $fq\n";
	system "gzip $fq\n";
	$pm->finish;
}

$pm->wait_all_children;



