#!/usr/bin/perl -w
use strict;
use Getopt::Std;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -f <RepeatMasker.out (additional .out files separated by commas)> > RepeatMasker.summary\n
Optional: -s Summary repeat name (r) or class (c) [c]
          -e Exclude low scoring matches within other matches (y/n) [y]\n";
our($opt_e, $opt_f, $opt_s);
getopt('efs');
die $usage unless ($opt_f);
if(!defined $opt_e) { $opt_e = 'y'; }
if(!defined $opt_s) { $opt_s = 'c'; }

# Summrise RepeatMasker
my %summary;

my @files = split /,/, $opt_f;
foreach my $file(@files) {
	warn "Saving summary for $file...\n";
	open my $fh, '<', $file or die "Cannot open $file: $!\n";
	while(my $line=<$fh>) {
		chomp $line;
		my $repeat_line = &parse_repeat_masker_out($line);

		# ignore blank, headers and low scoring matches within other matches
		next if($$repeat_line{'header'} eq 'Y');
		next if($$repeat_line{'next'} eq 1);
		if($opt_e ne 'n') {
			next if(defined $$repeat_line{'better_match_found'});
		}

		# Tally
		my $repeat_name = $$repeat_line{'repeat_name'};
		my $repeat_class = $$repeat_line{'repeat_class'};
		my $length = $$repeat_line{'length'};
		$summary{'repeat_name_nt'}{$repeat_name}{$file} += $length;
		$summary{'repeat_class_nt'}{$repeat_class}{$file} += $length;
		$summary{'total_nt'}{$file} += $length;
	}
	close $fh;
}

# Out summary header
my $type = 'repeat_name_nt';
if($opt_s eq 'r') { print "repeat_name"; }
else { 
	print "repeat_class"; 
	$type = 'repeat_class_nt';
}
foreach my $file(@files) { print "\t$file"; }
print "\n";

# Out summary
foreach my $repeat(sort keys %{$summary{$type}}) {
	print "$repeat";
	foreach my $file(@files) { 
		my $length = 0;
		if(defined $summary{$type}{$repeat}{$file}) { $length = $summary{$type}{$repeat}{$file}; }
		print "\t$length";
	}
	print "\n";
}

sub parse_repeat_masker_out {
	my $line = $_[0];
	my %line_info;

	# Ignore blank lines
	$line_info{'next'} = 0;
	$line_info{'header'}='N';
	if(($line =~ m/^$/) || ($line =~ m/^\n/)) {
		$line_info{'next'} = 1;
		return \%line_info;
	}
	# Ignore headers
	if(($line =~ m/\s+SW\s+perc\s+perc\s+perc/) || ($line =~ m/^score|\s+score/)) {
		$line_info{'next'} = 1;
		$line_info{'header'}='Y';
		return \%line_info;
	} 

	# Replace unusual Repeatmasker.out format
	$line =~ s/\s+/\t/g;
	if($line =~ m/^\t/) { $line =~ s/\t//;}
	my @bits = split /\t/, $line;
	$bits[7] =~ s/\(|\)//g;
	$bits[11] =~ s/\(|\)//g;
	$line_info{'smithwaterman_score'}	= $bits[0];
	$line_info{'percent_substitutions'}	= $bits[1];
	$line_info{'percent_deletions'}		= $bits[2];
	$line_info{'percent_insertions'}	= $bits[3];
	$line_info{'query_id'}			= $bits[4];
	$line_info{'start'}			= $bits[5];
	$line_info{'end'}			= $bits[6]; 
	$line_info{'bases_in_query_past_match'}	= $bits[7]; # no. of bases in query sequence past the ending position of match
	$line_info{'complement_in_db'}		= $bits[8];
	$line_info{'repeat_name'}		= $bits[9];
	$line_info{'repeat_class'}		= $bits[10];
	$line_info{'bases_in_repeat_past_match'}= $bits[11]; # no. of bases in (complement of) the repeat consensus sequence prior to beginning of the match (so 0 means that the match extended all the way to the end of the repeat consensus sequence)
	$line_info{'start_pos_in_db'}		= $bits[12]; # starting position of match in database sequence (using top-strand numbering)
	$line_info{'end_pos_in_db'}		= $bits[13]; # ending position of match in database sequence
	if(defined $bits[14]) { $line_info{'id_count'} = $bits[14]; }
	if(defined $bits[15]) {
		# An asterisk (*) in the final column (no example shown) indicates that there is a higher-scoring match whose domain partly (<80%) includes the domain of this match. 
		if($bits[15] eq '*') { $line_info{'better_match_found'} = 'Y'; }
	}
	die "6th column doesn't look right: $line" if($line_info{'end'} !~ m/\d+/);
	$line_info{'length'} = ($line_info{'end'} - $line_info{'start'});

	return \%line_info;
}
