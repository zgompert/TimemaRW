#!/usr/bin/env perl

# (c) Victor Soria-Carrasco
# victor.soria.carrasco@gmail.com
# Last modified: 29/10/2017 02:18:47

# Description:
# Converts between multiple alignment
# formats: fasta, nexus, and phylip
# (indicated with extension)

# Usage:
# aliconverter.pl -i <file> -o <file>

# Changelog
# 1.1a 
#	- Added option to output multistate alignment from DNA, RNA, or protein 
#	- Added option to output SNP (0,1,2) format from DNA (only if sites are biallelic) 
#	- Added automatic guessing of data type for nexus (dna, rna, protein, standard)
# 1.1b - 29/10/2017
#   - Output sequences sorted alphanumerically by id


use warnings;
use strict;
use File::Basename;
use File::Spec qw(rel2abs);
use Getopt::Long;
use Text::Wrap;
$Text::Wrap::columns=80;
use File::LocalizeNewlines; # handling of newlines for dos, mac, and unix

my $version='1.1b-2017.10.29';

&author;

my $ifile;
my $ofile;
my $transform='na';
my $datatype;
GetOptions(
    'i|I=s'       => \$ifile,
    'o|O=s'       => \$ofile,
    't|T=s'       => \$transform,
    'd|D=s'       => \$datatype,
    'h|help'      => \&usage
);

&usage if ((!defined $ifile) || (!defined $ofile));

my %formats=(
	'nex'=>'nexus',
	'nexus'=>'nexus',
	'phy'=>'phylip',
	'phylip'=>'phylip',
	'fa'=>'fasta',
	'fasta'=>'fasta',
	'fas'=>'fasta',
	'fal'=>'fasta',
	'fastal'=>'fasta');

$ifile=File::Spec->rel2abs($ifile);
my ($iname,$idir,$iext)=fileparse($ifile, keys %formats);
$ofile=File::Spec->rel2abs($ofile);
my ($oname,$odir,$oext)=fileparse($ofile, keys %formats);

&usage if (!defined($formats{$iext}) || !defined($formats{$oext}));

File::LocalizeNewlines->localize( $ifile );

my %ali;
print "Reading $formats{$iext} file $ifile...\n";
if ($formats{$iext} eq 'fasta'){
	%ali=read_fasta($ifile);
}
elsif ($formats{$iext} eq 'phylip'){
	%ali=read_phylip($ifile);
}
elsif ($formats{$iext} eq 'nexus'){
	%ali=read_nexus($ifile);
}

# Check it is an alignment
my $len=0;
foreach my $id (keys %ali){
	# print "$id - $len - ".length($ali{$id})." - $ali{$id}\n"; # debug
	die "\nThis is not a multiple alignment: sequences have different lenghts\n\n"
		if ($len !=0 && $len!=length($ali{$id}));
	$len=length($ali{$id});
}

# Convert it to multistate or SNP data if requested
if ($transform eq "m"){# Convert to multistate as in RAxML
	print "Converting to multistate data...\n";
	# up to 32 states as in RAxML
	
	%ali=convert2multistate(%ali);
}
elsif ($transform eq "s"){ # convert to SNP data (0,1,2) for SNAPP
	print "Converting to SNP data (0,1,2)\n";
	# 0-2 - counts of alternate allele
	# arbitrarily assign major allele as reference
	# ambiguous sites are interpreted as heterozygotes

	%ali=convert2snp(%ali);
}

my $ntaxa=scalar(keys %ali);
my $nsites=$len;

print "Converting to $formats{$oext} format...\n";
if ($formats{$oext} eq 'fasta'){
	write_fasta(\%ali, \$ofile);
	print "\nConverted file in fasta format saved as $ofile\n\n";
}
elsif ($formats{$oext} eq 'phylip'){
	write_phylip(\%ali, \$ofile, \$nsites);
	print "\nConverted file in phylip format saved as $ofile\n\n";
}
elsif ($formats{$oext} eq 'nexus'){
	write_nexus(\%ali, \$ofile, \$nsites);
	print "\nConverted file in nexus format saved as $ofile\n\n";
}


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
    print "  ".basename($0)." version $version     \n";
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
	print "      -t <convert to: multistate (m)| SNP (s) (optional)>\n";
	print "      -d <datatype: dna, rna, protein, standard (optional, try to guess otherwise)>\n";
	print "      -h <show this help>\n";
	print "";
	print "    Format is specified by extension (fasta,nexus,phy)\n";
    print "\n";
    exit;
}
# ==============================================================================

# Read fasta
# ==============================================================================
sub read_fasta{
	my $ifile=shift;

	open (FILE, "$ifile")
		or die "\n\nCan't open file $ifile\n\n";
		my $id='';
		my %ali=();
		while (<FILE>){
			chomp;
			if (/^\>/){
				s/^\>//g;
				$id=$_;
			}
			else{
				s/\s//g;
				$ali{$id}.=$_;
			}
		}
	close (FILE);

	return (%ali);
}
# ==============================================================================

# Read phylip
# ==============================================================================
sub read_phylip{ # requires separator between ids and sequences
	my $ifile=shift;

	open (FILE, "$ifile")
		or die "\n\nCan't open file $ifile\n\n";
		my $id='';
		my @ids=();
		my %ali=();
		my $block=0;
		my $c=0;
		while (<FILE>){
			# next if (/^[\t| ]*\n/);
			if (/^\s*\n/ && defined($ids[0])){
				$block++;
				$c=0;
			}
			elsif (!/[0-9]+\s+[0-9]+/){
				chomp;
				my @aux=split(/\s/,$_);
				my $seq='';
				if ($block>0){
					$id=$ids[$c];
					$seq.=$aux[0];
					$c++;
				}
				else{
					$id=$aux[0];
					push (@ids, $id);
				}
				$seq.=join('',@aux[1..$#aux]);
				$ali{$id}.=$seq;
			}
		}
	close (FILE);

	return (%ali);
}
# ==============================================================================

# Read nexus
# ==============================================================================
sub read_nexus{
	my $ifile=shift;

	open (FILE, "$ifile")
		or die "\n\nCan't open file $ifile\n\n";
		my $id='';
		my %ali=();
		my $read=0;
		while (<FILE>){
			if (/matrix/i){
				$read=1;
				next;
			}
			$read=0 if (/;/);
			next if (/^\[/ || /^[\s]*\n/);
			if ($read==1){
				chomp;
				s/\[.*//g;
				my @aux=split(/\s/,$_);
				$ali{$aux[0]}.=join('',@aux[1..$#aux]);
			}
		}
	close (FILE);

	return (%ali);
}
# ==============================================================================

# Write fasta
# ==============================================================================
sub write_fasta{
	my ($a, $of)=@_;
	my %ali=%$a;
	my $ifile=$$of;

	open (FILE, ">$ofile")
		or die "\n\nCan't write file $ofile\n\n";
		foreach my $id (sort keys %ali){
			print FILE ">$id\n";
			print FILE wrap ('','',$ali{$id}."\n");
		}
	close (FILE);
}
# ==============================================================================

# Write phylip
# ==============================================================================
sub write_phylip{
	my ($a, $of, $n)=@_;
	my %ali=%$a;
	my $ifile=$$of;
	my $nsites=$$n;
	
	my $ntaxa=scalar(keys %ali);
	open (FILE, ">$ofile")
		or die "\n\nCan't write file $ofile\n\n";
		print FILE " $ntaxa $nsites\n";
		foreach my $id (sort keys %ali){
			print FILE "$id\t$ali{$id}\n";
		}
	close (FILE);
}
# ==============================================================================

# Write nexus
# ==============================================================================
sub write_nexus{
	my ($a, $of, $n)=@_;
	my %ali=%$a;
	my $ifile=$$of;
	my $nsites=$$n;

	$datatype=guess_datatype(%ali)
		 if (!defined ($datatype));
	
	my $ntaxa=scalar(keys %ali);
	my %states;
	foreach my $id (keys %ali){
		my @aux=split('',$ali{$id});
		foreach my $a (@aux){ $states{$a}=1 if ($a ne '-' && $a ne '?'); }
	}
	my $symbols=join('',sort keys %states);

	open (FILE, ">$ofile")
		or die "\n\nCan't write file $ofile\n\n";
		print FILE "#NEXUS\n";
		print FILE "Begin data;\n";
		print FILE "Dimensions ntax=$ntaxa nchar=$nsites;\n";
		print FILE "Format datatype=$datatype gap=- missing=? symbols=\"$symbols\";\n";
		print FILE "Matrix\n";
		foreach my $id (sort keys %ali){
			print FILE "$id\t$ali{$id}\n";
		}
		print FILE ";\n";
		print FILE "End;\n";

	close (FILE);
}
# ==============================================================================

# Convert to multistate
# ==============================================================================
sub convert2multistate{
	my %ali=@_;

	my %states;
	# Get list of unique states
	foreach my $id (keys %ali){
		my %statesid=unique($ali{$id}); # get unique states for this sequence
		foreach my $s (keys %statesid){ $states{$s}+=$statesid{$s};	}

		my @aux=split(//,$ali{$id});
		foreach my $a (@aux){ $states{uc($a)}++ if ($a ne '-' || $a ne '?'); }
	}

	# Get equivalence to states (RAxML style)
	my @values=("0","1","2","3","4","5","6","7","8","9",
				"A","B","C","D","E","F","G","H","I","J",
				"K","L","M","N","O","P","Q","R","S","T",
				"U","V");
	
	die "\nERROR: Cannot encode as multistate, unique values in alignment are > 32\n\n"
		if (scalar (keys %states) > scalar(@values));

	# Assign states arbitrarily	
	my $i=0;
	foreach my $s (sort keys %states){
		$states{$s}=$values[$i];
		$i++;
	}

	# Replace states
	foreach my $id (keys %ali){
		foreach my $s (keys %states){
			$ali{$id}=~ s/$s/$states{$s}/gi;
		}
	}

	return(%ali);
}
# ==============================================================================

# Convert to SNP
# ==============================================================================
sub convert2snp{
	my %ali=@_;;
	# 0-2 - counts of alternate allele
	# arbitrarily assign major allele as reference
	# ambiguous sites are interpreted as heterozygotes
	
	# Get array of sites
	my @sites;
	my @ids;
	foreach my $id (sort keys %ali){
		my @sitesid=split(//,$ali{$id});
		foreach my $i (0..$#sitesid){
			$sites[$i].=$sitesid[$i];
		}
		push(@ids,$id);
	}
	%ali=();#empty alignment hash

	# Get unique states for each site, identify major allele, and replace with 0,1,2 values
	foreach my $i (0..$#sites){
		my %statessite=unique($sites[$i]);
		# check it is a trully biallelic SNP
		die "\nERROR: Sequence for pos $i contain more than 2 alleles (values: ".join(",",keys %statessite).")\n\n"
			if (scalar (keys %statessite) > 3 && # more than 3 states
				grep (/^[BDHVN]$/, keys %statessite)); # ambiguities involving 3+ alleles

		my $majal=0;
		foreach my $s (sort {$statessite{$b}<=>$statessite{$a}} keys %statessite){
			# print STDERR "$i - $s - $statessite{$s}\n";#debug
			if ($s=~ /[ACTG]/i){
				if ($majal==0){#homozygous for major/reference
					$statessite{$s}=0;
					$majal=1;
				}
				else{#homozygous for minor/alternate
					$statessite{$s}=2;
				}
			}
			else{#heterozygous
				$statessite{$s}=1;
			}
		}
		# Replace states
		foreach my $state (keys %statessite){
			$sites[$i]=~ s/$state/$statessite{$state}/ig;
		}
		my @statesid=split(//,$sites[$i]); # states for samples for this position
		foreach my $j (0..$#statesid){
			$ali{$ids[$j]}.=$statesid[$j];
		}
	}

	return(%ali);
}
# ==============================================================================

# Guess type of data
# ==============================================================================
sub guess_datatype{
	my %ali=@_;
	my $datatype='';
	my $dna='ACTGRYSWKMBDHVN';
	my $rna='ACUGRYSWKMBDHVN';
	my $protein='ACDEFGHIKLMNPQRSTVWYX';

	my %states;
	foreach my $id (keys %ali){
		my %statesid=unique($ali{$id}); # get unique states for this sequence
		foreach my $s (keys %statesid){ $states{$s}+=$statesid{$s};	}
	}
	
	my $isdna=my $isrna=my $isprotein=join('',keys %states);
	$isdna=~ s/[\-\?$dna]//ig;
	$isrna=~ s/[\-\?$rna]//ig;
	$isprotein=~ s/[\-\?$protein]//ig;
	
	my $check=0;
	if ($isdna eq ''){
		# print STDERR "DNA\n"; #debug;
		$datatype='dna';
		$check++;
	}
	if ($isrna eq ''){
		# print STDERR "RNA\n"; #debug;
		$datatype='rna';
		$check++;
	}
	if ($isprotein eq ''){
		# print STDERR "protein\n"; #debug;
		if ($check>0){ # putative DNA or RNA
			# check frequency ACTG
			# print STDERR "STATES: ".join('',keys %states)."\n";#debug
			my $total=0;
			foreach my $s (keys %states){$total+=$states{$s};}
			my $actgfreq=sprintf("%.2f",($states{'A'}+$states{'C'}+$states{'T'}+$states{'G'})/$total);
			# print STDERR "ACTGFR: $actgsum - $total - $actgfreq\n";#debug
			if ($actgfreq < 0.25){ # 25% is approx. expected frequency in vertebrates; it is lower in invertebrates
				$datatype='protein';
				print "\nWARNING: Couldn't guess type of data unambiguously, datatype set to 'protein'\n\n"
			}
			else{
				print "\nWARNING: Couldn't guess type of data unambiguously, datatype set to 'dna'\n\n"
			}
		}
		else{
			$datatype='protein';
			$check++;
		}
	}
	if ($check>1 || $check==0){
		$datatype='standard';
		print "\nWARNING: Couldn't guess type of data unambiguously, datatype set to 'standard'\n\n"
			if ($check>1);
	}
	return($datatype);
}
# ==============================================================================

# Get unique elements in string (and counts)
# ==============================================================================
sub unique{
	my $seq=shift;
	my @aux=split(//,$seq);
	my %elements;
	foreach my $a (@aux){
		$elements{uc($a)}++ 
			if ($a ne '-' && $a ne '?'); 
	}
	return (%elements);
}
# ==============================================================================

