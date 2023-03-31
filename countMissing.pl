#!/usr/bin/perl
#
open(IN, "perform_comb_og.fa");

while(<IN>){
	chomp;
	if(m/^>(\S+)/){
		$id = $1;
		$cnts{$id} = 0;
	}
	else{
		$c = () = $_ =~ /-/gi;
		$cnts{$id} = $c + $cnts{$id};
	}
}
close(IN);
foreach $id (sort keys %cnts){
	print "$id $cnts{$id}\n";
}
