#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use File::Basename;
#use FindBin qw($Bin);
#use lib "$Bin/perl_modules";

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -g <distances.tab from GFF-upstream-downstream-distances.pl>\n
Optional: -b\tNumber of bootstraps [5]
          -e\tEnrichment script [Density_plot_3and5prime_HgT.R]
	  -j\tJoin all stats script [Join_all_stats_bootstraps.pl]\n";
our($opt_b, $opt_e, $opt_g, $opt_j);
getopt('begj');
die $usage unless ($opt_g);
die "Cannot open $opt_g : $!" unless (-e $opt_g);
if(!defined $opt_b) { $opt_b = 5; }
if(!defined $opt_e) { $opt_e = 'Density_plot_3and5prime_HgT.R'; }
if(!defined $opt_j) { $opt_j = 'Join_all_stats_bootstraps.pl'; }
foreach($opt_e, $opt_j) { die "Can't find script : $_" if(! -e $_); }

# save file to memory
my ($distance_file_hash, $number_of_chosen_genes) = &save_distance_file($opt_g);
my $number_of_lines = scalar(keys(%{$distance_file_hash}));

# header
my $header = "contig\tfeature\tdesc\tstrand\tupstream_distance\tdownstream_distance\tstart_pos\tuplog10\tdownlog10";

# make new directory
my $dirname  = dirname($opt_g);
my $filename = fileparse($opt_g);
#die "dirname = $dirname\n";
my $outfolder = "$dirname/distance_bootstraps";
`mkdir $outfolder`;

# make random versions (cat1,) and test for enrichment for each
for(my $i=0; $i < $opt_b; $i++) {

	# Outfile
	my $outfile = ("$outfolder/$filename-random-bootstrapped-" . $i . '.tab');
	open my $ofh, '>', $outfile or die "Cannot open $outfile : $!";

	# Make random numbers
	my $random_numbers = &make_random_numbers_hash($number_of_lines, $number_of_chosen_genes);

	# Print header
	print $ofh "$header\n";

	# Go through file changing specified lines
	for(my $j=1; $j <= $number_of_lines; $j++) {
		my $line = $$distance_file_hash{$j};
		my $new_line = $line;

		# random number?
		if(defined $$random_numbers{$j}) {
			my @bits = split /\t/, $line;
			$bits[2] = 'cat1,';
			$new_line = join "\t", @bits;

		}

		# print
		print $ofh "$new_line\n";
	}
	close $ofh;

	# Density_plot_3and5prime_new_for_2speed_work_HgT
	my $cmd = "Rscript $opt_e -d $outfile";
	`$cmd`;
}

# join up stats for a summary file
warn "Join up stats...\n";
my $join_up_command = "perl $opt_j -g $outfolder";
`$join_up_command`;

# Outfiles
my $infile1 = "$outfolder/../stats_summary-HgT_Q_UL_p-values-bootstrap_random.txt";
my $infile2 = "$outfolder/../stats_summary-HgT_Q_UR_p-values-bootstrap_random.txt";
my $infile3 = "$outfolder/../stats_summary-HgT_Q_LR_p-values-bootstrap_random.txt";
my $infile4 = "$outfolder/../stats_summary-HgT_Q_LL_p-values-bootstrap_random.txt";

my $infile5 = "$outfolder/../stats_summary-X2_Q_UL_p-values-bootstrap_random.txt";
my $infile6 = "$outfolder/../stats_summary-X2_Q_UR_p-values-bootstrap_random.txt";
my $infile7 = "$outfolder/../stats_summary-X2_Q_LR_p-values-bootstrap_random.txt";
my $infile8 = "$outfolder/../stats_summary-X2_Q_LL_p-values-bootstrap_random.txt";

# summary outfile
warn "Generating summary file...\n";
my $summary_outfile = "$outfolder/../stats_summary-all.tab";
open my $ofh2, '>', $summary_outfile or die "Cannot open $summary_outfile : $!";

# Check for p < 0.01
foreach my $file($infile1, $infile2, $infile3, $infile4, $infile5, $infile6, $infile7, $infile8) {
	my $count_lines = 0;
	my $count_significance = 0;
	
	# open file
	open my $fh, '<', $file or die "Cannot open $file : $!";
	while(my $line=<$fh>) {
		chomp $line;
		my @bits = split /\t/, $line;
		if($bits[1] < 0.01) { $count_significance++; }
		$count_lines++;
	}

	# print report
	print $ofh2 "$file\t$count_lines\t$count_significance\n";
}
close $ofh2;

sub make_random_numbers_hash {
	my ($num_genes, $num_special_genes) = @_;
	my %random_numbers;

	while(scalar(keys(%random_numbers)) < $num_special_genes) {
		my $random_number = int(rand($num_genes));
		$random_numbers{$random_number} = 1;
	}
	my $sum_random_numbers = scalar(keys(%random_numbers));
	warn "make_random_numbers_hash: $sum_random_numbers made\n";

	return \%random_numbers;
}

sub save_distance_file {
	my $file = $_[0];
	my %file_hash;

	my $line_count = 1;
	my $number_of_cat1 = 0;
	open my $fh, '<', $file or die "Cannot open $file : $!";
	while(my $line=<$fh>) {
		chomp $line;

		# ignore header
		next if($line =~ m/^contig\tfeature/);

		# Save number of cat1 found, and then reset desc
		my @bits = split /\t/, $line;
		if($bits[2] eq 'cat1,') { $number_of_cat1++; }
		$bits[2] = 'NA';
		my $new_line = join "\t", @bits;

		# save
		$file_hash{$line_count} = $new_line;

		$line_count++;
	}
	warn "save_distance_file: $number_of_cat1 secreted / $line_count total found\n";
	return (\%file_hash, $number_of_cat1);
}
