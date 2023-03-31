#!/usr/bin/perl
#
# control script to estimate 2d SFS  with angsd
# use GL1 = samtools gl model 

## see
## http://popgen.dk/angsd/index.php/2d_SFS_Estimation


use Parallel::ForkManager;
my $max = 4;
my $pm = Parallel::ForkManager->new($max);

my $gl = 1;
my $genome = "/uufs/chpc.utah.edu/common/home/u6000989/data/timema/hic_genomes/t_knulli/mod_hic_output.fasta";


#-P 25 = 25 threads, with 3 jobs that is 75 threads

foreach $pop (@ARGV){
	$pm->start and next; ## fork
	$pop =~ m/^([a-zA-Z0-9\-_]+)\.fileslist/ or print "no match for $set\n";
	$base = $1;
	$o = "o_1d_$base";
	system "~/source/angsd/angsd -bam $pop -doSaf 1 -GL $gl -rf ScFnSIn_500_HRSCAF958:13093370-43606674 -P 25 -minMapQ 30 -minQ 20 -anc $genome -out $o\n";
	$pm->finish;
}

$pm->wait_all_children;

system "~/source/angsd/misc/realSFS $f1 $f2 -fold 1 -P 60 > $o1\n";

## now 2d
$f1 = "bce_cc.saf.idx";
$f2 = "bce_rw.saf.idx";
$f3 = "bcturn_cc.saf.idx";
system "~/source/angsd/misc/realSFS $f1 $f2 -fold 1 -P 60 > o2d_bce_cc_bce_rw\n";
system "~/source/angsd/misc/realSFS $f1 $f3 -fold 1 -P 60 > o2d_bce_cc_bcturn_cc\n";
system "~/source/angsd/misc/realSFS $f2 $f3 -fold 1 -P 60 > o2d_bce_rw_bcturn_cc\n";


