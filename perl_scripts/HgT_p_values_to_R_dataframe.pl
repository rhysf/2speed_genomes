#!/usr/bin/perl -w
use strict;
use Getopt::Std;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -a <HgT p-values QLL> -b <HgT p-values QLR> -c <HgT p-values QUL> -d <HgT p-values QUR> > outfile_dataframe.tab\n
Notes: Pay attention to make sure the quadrant name matches the name above. Input files are from Density_plot_3and5prime_new_for_2speed_work_HgT.R and Join_all_stats.pl\n";
our($opt_a, $opt_b, $opt_c, $opt_d);
getopt('abcd');
die $usage unless ($opt_a && $opt_b && $opt_c && $opt_d);
foreach($opt_a, $opt_b, $opt_c, $opt_d) { die "Cannot open $_ : $!" unless (-e $_); }

# Save values to arrays
my $QLL_array = &save_values($opt_a);
my $QLR_array = &save_values($opt_b);
my $QUL_array = &save_values($opt_c);
my $QUR_array = &save_values($opt_d);

# Print dataframe
print "Quadrant\tValue\n";
foreach my $value(@{$QLL_array}) {
	print "QLL\t$value\n";
}
foreach my $value(@{$QLR_array}) {
	print "QLR\t$value\n";
}
foreach my $value(@{$QUL_array}) {
	print "QUL\t$value\n";
}
foreach my $value(@{$QUR_array}) {
	print "QUR\t$value\n";
}

sub save_values {
	my $file = $_[0];
	my @values;

	open my $fh, '<', $file or die "Cannot open $file : $!";
	while(my $line=<$fh>) {
		chomp $line;
		next if($line =~ m/FILE\tp-value/);
		my @bits = split /\t/, $line;
		push @values, $bits[1];
	}
	close $fh;
	return \@values;
}