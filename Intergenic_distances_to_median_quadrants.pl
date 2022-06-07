#!/usr/bin/perl -w
use strict;
use Getopt::Std;
#use FindBin qw($Bin);
#use lib "$Bin/perl_modules";

### r.farrer@exeter.ac.uk

# Opening commands
my $usage = "Usage: perl $0 -i <output from GFF-upstream-downstream-distances.pl> > Quadrant_all_locations.tab\n";
our($opt_i);
getopt('i');
die $usage unless ($opt_i);
die "Cannot open $opt_i : $!" unless (-e $opt_i);

# Save uplog10 and downlog10 to array
my @uplog10 = &save_values_to_array($opt_i, 7);
my @downlog10 = &save_values_to_array($opt_i, 8);

# Calculate median
my $median_up = &median(@uplog10);
my $median_down = &median(@downlog10);
warn "median upstream (log10) = $median_up\n";
warn "median downstream (log10) = $median_down\n";

# Output genes and their quadrants
print "scaffold\tlocation\tquadrant\n";
open my $fh, '<', $opt_i or die "Cannot open $opt_i : $!\n";
while(my $line=<$fh>) {
	chomp $line;

	# ignore header
	next if($line=~ m/^contig\tfeature/);
	my @bits = split /\t/, $line;
	my ($contig, $feature, $desc, $strand, $upstream_distance, $downstream_distance, $start, $uplog10, $downlog10) = @bits;

	# quadrant?
	my $quadrant;
	if($downlog10 < $median_down) {
		if($uplog10 < $median_up) { $quadrant = 4; }
		else { $quadrant = 1; }
	} else {
		if($uplog10 < $median_up) { $quadrant = 3; }
		else { $quadrant = 2; }
	}
	die "ERROR: quadrant not defined for $uplog10 and $downlog10: $line\n" if(!defined $quadrant);

	# print
	print "$contig\t$start\t$quadrant\n";
}

sub median {
	my @vals = sort { $a <=> $b } @_;
	my $len = @vals;

	# odd
	if($len % 2) {
		return $vals[int($len / 2)];
	}
	# even
	else {
		return ($vals[int($len / 2) - 1] + $vals[int($len / 2)]) / 2;
	}
}

sub save_values_to_array {
	my ($file, $column) = @_;
	my @array;

	# open file and save column
	open my $fh, '<', $file or die "Cannot open $file : $!\n";
	while(my $line=<$fh>) {
		chomp $line;

		# ignore header
		next if($line=~ m/^contig\tfeature/);

		my @bits = split /\t/, $line;
		die "ERROR: Unexpected line (should have 8 columns from output from GFF-upstream-downstream-distances.pl): $line\n" if(!defined $bits[8]);
		die "ERROR: could not pull column $column from line $line\n" if(!defined $bits[$column]);
		push @array, $bits[$column];
	}
	return @array;
}
