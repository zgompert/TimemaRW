#!/usr/bin/perl
#
# subset to retain knulli CC and RW from BCE and petita
# also drop individuals with lots of missing data
# using 10% as cutoff = 79 missing sites of 789
# this saves 68% of individuals or 181

$tenpct = 79;

## get missingness counts
open(IN, "perform_comb_og.fa");
while(<IN>){
	chomp;
	if(m/^>t[a-z]+\-(\S+)/){
		$id = $1;
		$miss{$id} = 0;
	}
	else{
		$d = () = (m/-/g);
		$miss{$id} += $d;
	}
}
close(IN);

## now get list of BCE CC and BCE RW
## plus 8-10 randomly selected for each other species
foreach $file (@ARGV){
	$file =~ m/^([a-z_]+)\./;
	$group = $1;
	open(IN, $file) or die;
	while(<IN>){
		chomp;
		$gset{$_} = $group;
	}
	close(IN);
}		

## now subset alignment
open(IN, "perform_comb_og.fa");
open(OUT, "> sub_perform_comb_og.fa") or die;
$flag = 0;
while(<IN>){
	chomp;
	if(m/^>timema\-(\S+)/){
		$id = $1;
		if(($miss{$id} < $tenpct) and (exists($gset{$id}))){
			print "$id $miss{$id} $gset{$id}\n";
			$flag = 1;
			print OUT "$_\n";
		}
		else{
			$flag = 0;
		}
	}
	elsif(m/^>t[a-z]+\-(\S+)/){
		$id = $1;
		if(exists($gset{$id})){
			print "$id $miss{$id} $gset{$id}\n";
			$flag = 1;
			print OUT "$_\n";
		}
		else{
			$flag = 0;
		}
	}
	elsif($flag == 1){
		print OUT "$_\n";
	}
}
close(IN);
close(OUT);
