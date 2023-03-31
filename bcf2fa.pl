#!/usr/bin/env perl

# (c) Victor Soria-Carrasco
# victor.soria.carrasco@gmail.com
# Last modified: 15/02/2022 12:55:37

# Description:
# Converts a bcf file to a multiple aligment in fasta format
# See arguments below.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# Note: ignore positions with Ns in the reference
# ToDo:
# - Handle indels (now ignored)
# - Filter by likelihood value -> code as gap if genotype lnl < threshold
# - Compute consensus using a threshold

# Changelog
# 1.63 - 2017/11/28
#	- fixed bug on generating three-state output
#
# 1.62 - 2017/11/24
#	- fixed bug when using option 'strictvar'
#	- set to use specific version of bcftools
#	- enforced only snps are used (no indels)
#	- keep invariants if present
#	- removed bcftools option that prevented snps in reference Ns
#	- cleaned up code
#
# 1.61 - 2017/11/20
#	- fixed bug when using option 'strictvar'
#
# 1.6 - 2017/11/07
#   - randomised reference alleles for threestate output to avoid biases
#
# 1.51 - 2017/10/29
#   - fixed bug on generating three-state output
#
# 1.5 - 2017/05/11
#	- added option for three-state output (for SNAPP)
#	- some minor changes in on-screen help
#
# 1.4 - 2017/04/04
#   - added option 'strictvar' to exclude sites for which there are
#     not at least one homozygous sample for at least two different alleles
#
# 1.3 - 2017/02/24
#   - fixed serious bug that may cause loci to be written out in
#     different order for different samples
#
# 1.2 - 2016/08/02
#	- added option for multistate output

# ToDo - keep same order of loci (=chromosomes/scaffolds) as in the original file

use warnings;
use strict;
use File::Basename;
use File::Spec qw(rel2abs);
use Getopt::Long;
use List::MoreUtils qw(uniq);
use Text::Wrap;
$Text::Wrap::columns=80;

my $version='1.63-2017.11.28';

&author;

# External programs:
my $bcftools='bcftools';
# my $bcftools='/usr/local/extras/Genomics/apps/bcftools/1.5/bin/bcftools';
my $bcftools_version=`$bcftools |& grep Version | awk '{print \$2}'`;

my $ifile; # input file
my $ofile; # output file
my $rllist; # list to rename loci
my $rslist; # list to rename samples
my $maxpgaps='1.0'; # proportion of gaps per site allowed
my $maxphets='1.0'; # proportion of heterozyogotes per site allowed
my $maxpgaphets='1.0'; # proportion of heterozyogotes per site allowed

my $maxpgapt='1.0'; # proportion of gaps per taxa allowed
my $maxphett='1.0'; # proportion of heterozyogotes per taxa allowed
my $maxpgaphett='1.0'; # proportion of heterozyogotes per taxa allowed

my $strictvar=0; # require at least one homozygote for at least two different allele
my $genolhl=0; # use PLs to identify genotypes instead of GT
my $partition=0; # partition by locus
my $multistate=0; # output multistate instead
my $threestate=0; # output threestate instead (SNAPP)
my $consensus=0; # make consensus pooling multiple samples together
my $consfield=0; # field on samples names to pool samples together
my $conssplchr='-'; # character separating fields

my $refseq="/uufs/chpc.utah.edu/common/home/u6000989/data/timema/hic_genomes/t_knulli/mod_hic_output.fasta"; # reference sequence, required for 3-state SNAPP encoding

GetOptions(
    'i|I=s'       => \$ifile,
    'o|O=s'       => \$ofile,
	'p|P'         => \$partition,
	'm|M'         => \$multistate,
	't|T'         => \$threestate,
	'c|C'         => \$consensus,
	's|S'         => \$strictvar,
	'g|G'         => \$genolhl,
	'cf|CF=i'     => \$consfield,
	'rl|RL=s'     => \$rllist,
	'rs|RS=s'     => \$rslist,
	'pgs|PGS=f'   => \$maxpgaps,
	'phs|PHS=f'   => \$maxphets,
	'pghs|PGHS=f' => \$maxpgaphets,
    'pgt|PGT=f'   => \$maxpgapt,
	'pht|PHT=f'   => \$maxphett,
	'pght|PGHT=f' => \$maxpgaphett,
	'h|help'      => \&usage
);

&usage if	(!defined($ifile) || 
			($partition==1 && (!defined($ofile))) ||
			($consensus==1 && (!defined($rslist)))) ;

$ifile=File::Spec->rel2abs($ifile);
my $input="$ifile";
if ($input =~ m/\.bcf$/i){
	if ($bcftools_version=~ m/^0\.1/){ # bcftools 0.1.x
		# only snvs (-I to skip indels)
		# $input="$bcftools view -I -N $ifile |";
		$input="$bcftools view -I $ifile |";
	}
	else{ # bcftools v 1.x
		# only snvs (-v snps)
		# $input="$bcftools view -V indels -e 'REF=\"N\"' $ifile |";
		$input="$bcftools view -V indels $ifile |";
	}
}
elsif ($input =~ m/\.vcf/i){
	$input="cat $ifile |";
}
else{
	die ("\nFormat of input file $ifile not recognised\n\n");
}

# Read list of loci
my %rll;
if (defined($rllist)){
	open (FILE, "$rllist")
		or die ("\nCan't open file $rllist\n\n");
		while (<FILE>){
			chomp;
			my @aux=split(/\s|\,/,$_);
			$rll{$aux[0]}=$aux[1];
		}
	close (FILE);
}

# Read list of samples
my %rsl;
if (defined($rslist)){
	open (FILE, "$rslist")
		or die ("\nCan't open file $rslist\n\n");
		while (<FILE>){
			chomp;
			my @aux=split(/\s|\,/,$_);
			$rsl{$aux[0]}=$aux[1];
		}
	close (FILE);
}

my %seqs;
my %refseqs;
my %lociorder;
my @ids;
# my $nsnps=`$input grep -chv "^#"`;
# chomp($nsnps);
my $n=0;
my $c=0;
my $lo=0;
open (FILE, $input);
	while (<FILE>){
		if (/^\#/){ # header
			if (/CHROM/){
				chomp;
				my @aux=split(/\s/,$_);
				@ids=@aux[9..$#aux];
				if (defined($rslist)){
					foreach my $id (@ids){
						$id=~ s/(\.sorted)?(\.bam)?$//g;
						# print "$id -> $rsl{$id}\n";
						$id=$rsl{$id}."-".$id;
					}
				}
			}
		}
		else{
			$n++;
			chomp;
			my @aux=split(/\t/,$_);
			my @alleles;
			my $locus=$aux[0];
			$locus=$rll{$locus} if (defined($rllist));
			
			if (!exists($lociorder{$locus})){
				$lociorder{$locus}=$lo;
				$lo++;
			}

			$alleles[0]=$aux[3]; # reference allele
			push (@alleles, split(/\,/,$aux[4])); # alternate alleles
			
			my $refal=0;
			# random reference allele for three states encoding (to avoid biases)
			$refal=int(rand @alleles) if ($threestate==1);
			
			my $pgaps=0;
			my $phets=0;
			# my %sites=();#debug
			# my $seq=''; #debug
			if ($genolhl==0){ # do not use genotype likelihoods
				for my $i (9..$#aux){ # samples
					
					$refseqs{$ids[$i-9]}{$locus}.=$alleles[$refal];

					# gaps
					if (($aux[$i]=~ m/\:(0\,0\,0)+/) || # no data, coded as GL or PL = 0,0,0
						 $aux[$i]=~ m/^0$/ || # special case for bcftools 0.1x, no data coded as a single 0
						 $aux[$i]=~ m/^\.\/?\.?\:/){ # generic case (for bcftools 1.x, freebayes, etc)
							$seqs{$ids[$i-9]}{$locus}.='-';
							# $seq.='-';#debug
							$pgaps++;
					}
					elsif ($aux[$i] =~ /^([0-3])\|([0-3])\:/){ # phased
						$phets++ if ($alleles[$1] ne $alleles[$2]);
						
						if ($threestate==1 || $multistate==1){ # force use of IUPAC ambiguity codes
							if ($alleles[$1] eq $alleles[$2]){
								$seqs{$ids[$i-9]}{$locus}.=$alleles[$1];
							}
							else{
								$seqs{$ids[$i-9]}{$locus}.=code_heterozygote($alleles[$1].$alleles[$2]);
							}
						}
						else{ # get first haplotype
							$seqs{$ids[$i-9]}{$locus}.=$alleles[$1];
						}
					}
					elsif ($aux[$i] =~ /^([0-3])\/([0-3])\:/){ # unphased use IUPAC ambiguity codes
						$phets++ if ($alleles[$1] ne $alleles[$2]);

						if ($alleles[$1] eq $alleles[$2]){
							$seqs{$ids[$i-9]}{$locus}.=$alleles[$1];
						}
						else{
							$seqs{$ids[$i-9]}{$locus}.=code_heterozygote($alleles[$1].$alleles[$2]);
						}
					}
				}
			}
			else{
			}
			
			# Number of alleles with at least one homozygous sample
			my $nstates=0;
			if ($strictvar==1){
				my $seq="";
				for my $i (9..$#aux){ $seq.=substr($seqs{$ids[$i-9]}{$locus},-1);} # states for last added site
				my $uniquestates=join("",uniq(split(//,$seq)));
				$nstates=()=$uniquestates=~ /[actg]/gi;
			}

			if ($maxpgaps <1 || $maxphets < 1 || $maxpgaphets < 1){
				# Calculate combined proportion of gaps and heterozygotes for this site
				my $pgaphets=($pgaps+$phets)/(scalar(@aux)-9);	
				# Calculate proportion of gaps for this site
				$pgaps=$pgaps/(scalar(@aux)-9);
				# Calculate proportion of heterozygotes for this site
				$phets=$phets/(scalar(@aux)-9);


				if ($pgaps > $maxpgaps || $phets > $maxphets || $pgaphets > $maxpgaphets ||
					($strictvar==1 && $nstates < 2)){
					for my $i (9..$#aux){
						chop($seqs{$ids[$i-9]}{$locus});
						chop($refseqs{$ids[$i-9]}{$locus});
					}
				}
				else{
					# print "alleles: ".join("/", @alleles)."\n";#debug
					# print "Retaining site $n - $aux[0]:$aux[1] ".join("/",keys %sites)." - $seq\n";#debug
					$c++;
				}
			}
			elsif ($strictvar==1 && $nstates < 2){
				for my $i (9..$#aux){
					# chop($seqs{$ids[$i-9]});
					chop($seqs{$ids[$i-9]}{$locus});
					chop($refseqs{$ids[$i-9]}{$locus});
				}
			}
			else{
				# print "Retaining site $n - ".join("/",keys %sites)."\n";#debug
				$c++;
			}
			
			print STDERR " # SNPs/sites retained: $c out of $n\r" if ($n%1000==0 || eof);
			# $|++;
		}
	}
close (FILE);

print STDERR " # SNVs/sites retained: $c out of $n\r";
print STDERR "\n\n";



my @loci=sort keys %{$seqs{$ids[0]}};


# Do consensus
# ------------------------------------------------------------------------------
if ($consensus==1){
	# create arrays pooling samples
	my %taxa;
	foreach my $id (sort keys %seqs){
		my @aux=split(/$conssplchr/,$id);
		push(@{$taxa{$aux[$consfield]}},$id);
	}
	
	# compute consensus sequence for each taxon (= pooled seqs)
	my %newseqs;
	foreach my $t (sort keys %taxa){
		foreach my $l (@loci){
			my @seqs=();
			foreach my $id (@{$taxa{$t}}){
				push (@seqs, $seqs{$id}{$l});
				$|++;
			}
			$newseqs{$t}{$l}=consensus(@seqs);
		}
	}
	
	%seqs=%newseqs;
}
# ------------------------------------------------------------------------------6

# Output
# ------------------------------------------------------------------------------
my %locus_files;
if (defined($ofile)){
	if ($partition==0){
		open (FILEOUT, ">$ofile") 
			or die "(\nCan't write to output file $ofile\n\n)";
	}
	else{
		# open a file for each locus
		foreach my $l (@loci){
			my $of=$ofile;
			$of=~ s/\.fa(sta)?//g;
			$of.="-$l.fa";
			open (my $fh, "> $of")
				or die ("\nCan't write to $of\n\n");
			$locus_files{$l}=$fh;
		}
	}
}
else{
	*FILEOUT = *STDOUT;
}

if ($maxpgapt < 1 || $maxphett < 1 || $maxpgaphett < 1){
	foreach my $id (sort keys %seqs){
		# Calculate proportion of gaps and heterozygotes for each sample/taxon
		# Concatenate all loci
		my $seq='';
		my $refseq='';
		foreach my $l (sort keys %{$seqs{$id}}){
			$seq.=$seqs{$id}{$l};
			$refseq.=$refseqs{$id}{$l};
		}

		my $pgaphett=0;
		
		# Calculate proportion of gaps for this sample/taxon
		my $pgapt=()=$seq =~ /\-/g;
		$pgaphett+=$pgapt;
		$pgapt/=length($seq);

		# Calculate proportion of heterozygotes for this sample/taxon
		my $phett=()=$seq=~ /[KMRSWY]/ig;
		$pgaphett+=$phett;
		$phett/=length($seq);

		# Calculate combined proportion of gaps and heterozygotes for this sample/taxon
		$pgaphett/=length($seq);
		
		if ($pgapt <= $maxpgapt && $phett <= $maxphett && $pgaphett <= $maxpgaphett){
			my $outid=$id;
			$outid=~ s/(\.sorted)?\.bam//g; # remove extensions
			if ($partition==0){
				print FILEOUT ">$outid\n";
				$seq=code_multistate($seq) if ($multistate==1);
				$seq=code_threestate($refseq,$seq) if ($threestate==1);
				print FILEOUT wrap('','', $seq."\n");
			}
			else{
				#Output each locus to a separate file
				foreach my $l (sort keys %{$seqs{$id}}){
					print {$locus_files{$l}} ">$outid\n";
					$seqs{$id}{$l}=code_multistate($seqs{$id}{$l}) if ($multistate==1);
					$seqs{$id}{$l}=code_threestate($refseqs{$id}{$l}, $seqs{$id}{$l}) if ($threestate==1);
					print {$locus_files{$l}} wrap('','', $seqs{$id}{$l}."\n");
				}
			}
		}
		else{
			print STDERR "$id excluded (pgap: $pgapt > $maxpgapt or phet: $phett > $maxphett, or pgaphet: $pgaphett > $maxpgaphett)\n";
		}
	}
}
else{
	foreach my $id (sort keys %seqs){
		my $outid=$id;
		$outid=~ s/(\.sorted)?\.bam//g; # remove extensions

		if ($partition==0){
			# Concatenate all loci
			my $seq='';
			my $refseq='';
			foreach my $l (sort {$lociorder{$a}<=>$lociorder{$b}} keys %lociorder){
				# print "$id - $l - $seqs{$id}{$l} $refseqs{$id}{$l}\n"; # debug
				$seq.=$seqs{$id}{$l};
				$refseq.=$refseqs{$id}{$l};
			}
			print FILEOUT ">$outid\n";
			$seq=code_multistate($seq) if ($multistate==1);
			$seq=code_threestate($refseq, $seq) if ($threestate==1);
			print FILEOUT wrap('','', $seq."\n");
		}
		else{
			#Output each locus to a separate file
			foreach my $l (sort {$lociorder{$a}<=>$lociorder{$b}} keys %lociorder){
				print {$locus_files{$l}} ">$outid\n";
				$seqs{$id}{$l}=code_multistate($seqs{$id}{$l}) if ($multistate==1);
				$seqs{$id}{$l}=code_threestate($refseqs{$id}{$l}, $seqs{$id}{$l}) if ($threestate==1);
				print {$locus_files{$l}} wrap('','', $seqs{$id}{$l}."\n");
			}
		}
	}
}

if ($partition==0){
	close (FILEOUT) if (defined($ofile));
}
else{
	foreach my $l (@loci){
		close ($locus_files{$l});
	}
}
# ------------------------------------------------------------------------------

# ==============================================================================
# ==============================================================================
# ============================== SUBROUTINES ===================================
# ==============================================================================
# ==============================================================================


# Show copyright
# ==============================================================================
sub author{
    print "\n";
    print "#########################################\n";
    print "  ".basename($0)."\n";
	print "  version $version     \n";
    print "  (c) Victor Soria-Carrasco             \n";
    print "  victor.soria.carrasco\@gmail.com      \n";
    print "#########################################\n";
	print "\n";
}
# ==============================================================================

# Show usage
# ==============================================================================
sub usage{
    print "\n";
	print "  Usage:\n";
    print "    ".basename($0)."\n";
	print "      -i <input file>\n";
	print "      -o <output file>\n";
	print "      -p (output each locus in a separate file; requires specifying an output file)\n";
	print "      -m (output multistate (32 states: [0-9A-V]) encoding)\n";
	print "      -t (output threestate (0/1/2) SNAPP encoding)\n";
	print "      -c (make consensus pooling multiple samples; requires specifying a list of samples with -rs )\n";
	print "      -s (strict variant mode; output a variant site only if there is at least one homozygous sample for each allele)\n";
	print "      -rl   <csv file> (rename loci using this list)\n";
	print "      -rs   <csv file> (rename samples using this list)\n";
	print "      -pgs  <proportion of gaps per site allowed (optional; 0-1; default:1.0)>\n";
	print "      -phs  <proportion of heterozygotes per site allowed (optional; 0-1; default:1.0)>\n";
	print "      -pghs <combined proportion of gaps and heterozygotes per site allowed (optional; 0-1; default:1.0)>\n";
   	print "      -pgt  <proportion of gaps per taxa allowed (optional; 0-1; default:1.0)>\n";
	print "      -pht  <proportion of heterozygotes per taxa allowed (optional; 0-1; default:1.0)>\n";
	print "      -pght <combined proportion of gaps and heterozygotes per taxa allowed (optional; 0-1; default:1.0)>\n";
    print "\n";
    exit;
}
# ==============================================================================

# Code heterozygotes according to IUPAC
# ==============================================================================
sub code_heterozygote{
	my %het=(
		'AC'=>'M',
		'AG'=>'R',
		'AT'=>'W', 
		'CA'=>'M', 
		'CG'=>'S', 
		'CT'=>'Y',
		'GA'=>'R', 
		'GC'=>'S', 
		'GT'=>'K',
		'TA'=>'W',
		'TC'=>'Y',
		'TG'=>'K'
	);

	return ($het{$_[0]});
}
# ==============================================================================

# Resolve heterozygotes encodedaccording to IUPAC
# ==============================================================================
sub rev_code_heterozygote{
	my %revhet=(
		'M'=>'AC',
		'R'=>'AG',
		'W'=>'AT', 
		'S'=>'CG', 
		'Y'=>'CT',
		'K'=>'GT'
	);
	my @out=split('',$revhet{$_[0]});
	return (@out);
}
# ==============================================================================

# IUPAC ambiguities
# ==============================================================================
sub code_ambiguities{
	
	my %amb=(
		'AC'=>'M', 'CA'=>'M', 
		'AG'=>'R', 'GA'=>'R',
		'AT'=>'W', 'TA'=>'W',
		'CG'=>'S', 'GC'=>'S',
		'CT'=>'Y', 'TC'=>'Y',
		'GT'=>'K', 'TG'=>'K',
		'ACG'=>'V',	'AGC'=>'V',	'CGA'=>'V',	'CAG'=>'V',	'GAC'=>'V',	'GCA'=>'V',
		'ACT'=>'H',	'ATC'=>'H',	'CTA'=>'H',	'CAT'=>'H',	'TAC'=>'H',	'TCA'=>'H',
		'AGT'=>'D',	'ATG'=>'D',	'GTA'=>'D',	'GAT'=>'D',	'TAG'=>'D',	'TGA'=>'D',
		'CGT'=>'B',	'CTG'=>'B',	'GTC'=>'B',	'GCT'=>'B',	'TCG'=>'B',	'TGC'=>'B',
		'ACGT'=>'N', 'ACTG'=>'N', 'ATCG'=>'N', 'TACG'=>'N', 'TAGC'=>'N', 'ATGC'=>'N',
		'AGTC'=>'N', 'AGCT'=>'N', 'GACT'=>'N', 'GATC'=>'N', 'GTAC'=>'N', 'TGAC'=>'N',
		'TGCA'=>'N', 'GTCA'=>'N', 'GCTA'=>'N', 'GCAT'=>'N', 'CGAT'=>'N', 'CGTA'=>'N',
		'CTGA'=>'N', 'TCGA'=>'N', 'TCAG'=>'N', 'CTAG'=>'N', 'CATG'=>'N', 'CAGT'=>'N'
	);

	return ($amb{$_[0]});
}
# ==============================================================================

# Code threestate for SNAPP
# ==============================================================================
sub code_threestate{
	# 0 homozygote for reference 
	# 1 heterozygote
	# 2 homozygote for alternate
	my $refseq=shift;
	my $seq=shift;
	my @aux1=split(//,$refseq);
	my @aux2=split(//,$seq);
	my %states=(
		'A'=>'2',
		'C'=>'2',
		'T'=>'2',
		'G'=>'2',
		'M'=>'1',
		'R'=>'1',
		'W'=>'1',
		'S'=>'1',
		'Y'=>'1',
		'K'=>'1',
		'V'=>'?',
		'H'=>'?',
		'D'=>'?',
		'B'=>'?',
		'N'=>'?'
	);

	$seq='';
	foreach my $i (0..$#aux1){
		if ($aux2[$i] eq $aux1[$i]){
			$seq.=0;	
		}
		else{
			if (defined($states{$aux2[$i]})){
				$seq.=$states{$aux2[$i]};
			}
			else{
				$seq.='?';
			}
		}
	}
	
	return ($seq);
}
# ==============================================================================

# Code multistate for RAxML
# ==============================================================================
sub code_multistate{
	my $seq=shift;
	my @aux=split(//,$seq);
	my %state=(
		'A'=>'0',
		'C'=>'1',
		'T'=>'2',
		'G'=>'3',
		'M'=>'4',
		'R'=>'5',
		'W'=>'6',
		'S'=>'7',
		'Y'=>'8',
		'K'=>'9',
		'V'=>'A',
		'H'=>'B',
		'D'=>'C',
		'B'=>'D',
		'N'=>'E'
	);

	$seq='';
	foreach my $a (@aux){
		if (defined($state{$a})){
			$seq.=$state{$a};
		}
		else{
			$seq.=$a;
		}
	}
	
	return ($seq);
}
# ==============================================================================

# Do consensus
# ==============================================================================
sub consensus{
	my @ali=@_;
	my $conseq='';
	my $thr=0.5; # threshold to use ambiguities (NOT IMPLEMENTED)
	while (length($ali[0]) > 0){
		my %cnt;
		my $totcnt=0;
		foreach my $s (@ali){
			my $nt=substr ($s, 0, 1, "");
			$nt=uc($nt);
			if ($nt !~ /[A|T|C|G|-]/){
				my @nts=rev_code_heterozygote($nt);
				foreach my $n (@nts){
					$cnt{$n}++;
				}
			}
			else{
				$cnt{$nt}+=2;
			}
			$totcnt+=2;
		}
		my $max=0;
		my $fnt='-';
		foreach my $nt (keys %cnt){
			if ($nt ne '-' && $cnt{$nt} > $max){
				$fnt=$nt;
				$max=$cnt{$nt};
			}
			elsif ($nt ne '-' && $cnt{$nt} == $max){
				$fnt=code_heterozygote($fnt.$nt);
			}
		}
		$conseq.=$fnt;

	}
	return($conseq);
}
# ==============================================================================


