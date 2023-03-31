#!/usr/bin/perl
#
# splits a genotype likelihood file by population and host
#
# usage splitPops.pl in.gl
#


open(TAB, "pop_plant_location_ids.csv") or die "failed to read the ID table\n";
<TAB>; ## burn header
while(<TAB>){
	chomp;
	@line = split(",",$_);
	$pid{$line[2]} = "$line[4]_$line[3]";
	print "$line[2] = $line[4]_$line[3]\n";
}

close(TAB);

my $in = shift (@ARGV); ## gl file

$in =~ m/(t[a-zA-z]+)\.gl/ or die "no match\n";
$base = $1;

open(IN, $in) or die "failed to open the infile\n";

## get no. loci, no. inds
$line = <IN>;
chomp($line);
$line =~ m/(\d+)\s+(\d+)/;
$nind = $1;
$nloc = $2;

## get ind. and pop. ids
$line = <IN>;
chomp($line);
@line = split (" ",$line);
foreach $ind (@line){
	$ind =~ m/(SRR\d+)/;
	$id = $pid{$1};
	print "$ind -- $id\n";
	push (@id,$id);
	push (@{$popids{$id}},$ind);
	$ids{$id} = 1;
	if(defined $popn{$id}){
		$popn{$id}++;
	}
	else {
		$popn{$id} = 1;
	}
}

## open one file per population
foreach $id (sort keys %ids){
	$fh = "F"."$id";
	open ($fh, "> $base"."_$id".".gl") or die "Could not write $id\n";
	$files{$id} = $fh;
	print {$files{$id}} "$popn{$id} $nloc\n";
	$pids = join (" ",@{$popids{$id}});
	print {$files{$id}} "$pids\n";
	@ones = ();
	for($i=0;$i<$popn{$id}; $i++){
		push (@ones,1);
	}
	$ones = join (" ", @ones);
	print {$files{$id}} "$ones\n";
}

## read and write
while (<IN>){
	chomp;
	@line = split (" ",$_);
	$a = shift(@line); ## locus info
	foreach $id (sort keys %ids){
		print {$files{$id}} "$a";
	}
	for ($i=0; $i<$nind; $i++){
		$id = $id[$i];
		for ($j=0; $j<3; $j++){
			$a = shift(@line); 
			print {$files{$id}} " $a";
		}
	}	
	foreach $id (sort keys %ids){
		print {$files{$id}} "\n";
	}
			
}
close (IN);
