#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/perl_modules";
use Data::Dumper;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -d <download_fna_files.sh> -g <folder of GFFs downloaded from NCBI>\n
Notes: PRrobably need to use from same folder as -g\n";
our($opt_d, $opt_g);
getopt('dg');
die $usage unless ($opt_d && $opt_g);
die "Cannot open $opt_d : $!" unless (-e $opt_d);

# Save all GFF files
warn "Save each Gff file in $opt_g (suffix .gff)...\n";
my @gff_files = <$opt_g/*.gff>;

# Loop through all gff files, saving accessions, and making a directory of that accession to move the gff to.
my %accessions_saved;
GFF: foreach my $gff(@gff_files) {

	# Next look up TaxID (e.g. .//GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.gff)
	my @gff_parts = split /\_/, $gff;
	$gff_parts[0] =~ s/\.\/+//;
	my $accession = ($gff_parts[0] . '_' . $gff_parts[1]);
	$accessions_saved{$accession} = 1;

	# mkdir for them, and move gff to that directory!
	`mkdir $accession`;
	`mv $gff $accession/`;
}
warn "I've found this number of accessions/GFFs = " . scalar(keys(%accessions_saved)) . "\n";

# Go through download_fna_files.sh and download those genome files

my $matched = 0;
open my $fh, '<', $opt_d or die "CAnnot open $opt_d : $!\n";
while(my $line=<$fh>) {
	chomp $line;
	foreach my $accession(keys %accessions_saved) {
		next if($line !~ m/$accession/);

		# found a corresponding GFF
		$matched++;
		if($matched % 100) { }
		else { warn "So far found $matched\n"; }
		delete $accessions_saved{$accession};

		# download genome file and move to file
	
		# download
		system($line);
	
		# move it
		`mv *.gz $accession/`;
	}
}
warn "match = $matched\n";
warn "not matched = " . scalar(keys(%accessions_saved)) . "\n";
foreach my $accession(keys %accessions_saved) {
	print "$accession\n";
}
