#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/perl_modules";
use MutationTools::read_GFF;
use Data::Dumper;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -g <directory of assembly directories with 1 Enrichment_tests.tab file in each>\n
Notes: Creates 8 files in the directory behind opt_g\n";
our($opt_g);
getopt('g');
die $usage unless ($opt_g);

# outfiles
my $outfile1 = "$opt_g/../stats_summary-HgT_Q_UL_p-values-bootstrap_random.txt";
my $outfile2 = "$opt_g/../stats_summary-HgT_Q_UR_p-values-bootstrap_random.txt";
my $outfile3 = "$opt_g/../stats_summary-HgT_Q_LR_p-values-bootstrap_random.txt";
my $outfile4 = "$opt_g/../stats_summary-HgT_Q_LL_p-values-bootstrap_random.txt";

my $outfile5 = "$opt_g/../stats_summary-X2_Q_UL_p-values-bootstrap_random.txt";
my $outfile6 = "$opt_g/../stats_summary-X2_Q_UR_p-values-bootstrap_random.txt";
my $outfile7 = "$opt_g/../stats_summary-X2_Q_LR_p-values-bootstrap_random.txt";
my $outfile8 = "$opt_g/../stats_summary-X2_Q_LL_p-values-bootstrap_random.txt";

open my $ofh1, '>', $outfile1 or die "Cannot open $outfile1";
open my $ofh2, '>', $outfile2 or die "Cannot open $outfile2";
open my $ofh3, '>', $outfile3 or die "Cannot open $outfile3";
open my $ofh4, '>', $outfile4 or die "Cannot open $outfile4";

open my $ofh5, '>', $outfile5 or die "Cannot open $outfile5";
open my $ofh6, '>', $outfile6 or die "Cannot open $outfile6";
open my $ofh7, '>', $outfile7 or die "Cannot open $outfile7";
open my $ofh8, '>', $outfile8 or die "Cannot open $outfile8";

# Save all distance files (changed from $opt_g/*/*Enrichment_tests.tab which was used for all genomes. Changed for bootstrap analysis)
warn "Save each stats file in $opt_g (suffix distances.tab)...\n";
#my @enrichment_test_files = <$opt_g/*/*Enrichment_tests.tab>;
my @enrichment_test_files = <$opt_g/*Enrichment_tests.tab>;


# Loop through all gff files
FILES: foreach my $enrich_file(@enrichment_test_files) {


	# print and save contig names
	warn "opening $enrich_file...\n";
	open my $fh, '<', $enrich_file or die "Cannot open $enrich_file : $!";
	while(my $line=<$fh>) {
		chomp $line;

		# ignore blank lines
		next if($line =~ m/^\n/);
		next if($line =~ m/^$/);

		my @bits = split /\t/, $line;
		my ($test, $value) = @bits;

		# save
		if($test eq 'cat1 p(hypergeometric) Q_UL:') { print $ofh1 "$enrich_file\t$value\n"; }
		if($test eq 'cat1 p(hypergeometric) Q_UR:') { print $ofh2 "$enrich_file\t$value\n"; }
		if($test eq 'cat1 p(hypergeometric) Q_LR:') { print $ofh3 "$enrich_file\t$value\n"; }
		if($test eq 'cat1 p(hypergeometric) Q_LL:') { print $ofh4 "$enrich_file\t$value\n"; }

		if($test eq 'cat1 p(chi-squared) Q_UL:') { print $ofh5 "$enrich_file\t$value\n"; }
		if($test eq 'cat1 p(chi-squared) Q_UR:') { print $ofh6 "$enrich_file\t$value\n"; }
		if($test eq 'cat1 p(chi-squared) Q_LR:') { print $ofh7 "$enrich_file\t$value\n"; }
		if($test eq 'cat1 p(chi-squared) Q_LL:') { print $ofh8 "$enrich_file\t$value\n"; }
	}
	#die "end here!\n";
}
close $ofh1;
close $ofh2;
close $ofh3;
close $ofh4;
close $ofh5;
close $ofh6;
close $ofh7;
close $ofh8;
