#!/usr/bin/perl
#
# subset only the perform locus
# ScFnSIn_500_HRSCAF958:13093370-43606674
#
#

$vcf = shift(@ARGV);
$cnt = 0;
open(IN, $vcf) or die;
open(OUT, "> perform_comb.vcf") or die;
while(<IN>){
	chomp;
	if(m/^Sc/){
		if(m/^ScFnSIn_500_HRSCAF958\s+(\d+)/){
			$pos = $1;
			if(($pos >= 13093370) and ($pos <= 43606674)){
				$cnt++;
				print OUT "$_\n";
			}
		}
	}
	else{
		print OUT "$_\n";
	}
}
close(IN);
close(OUT);
print "Retained $cnt SNPs\n";
