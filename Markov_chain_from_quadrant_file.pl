#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Math::Matrix;
use Data::Dumper;

### r.farrer@exeter.ac.uk

# usage
my $usage = "perl $0 -q <quadrant_all_locations.tab (scaffold location quadrant)>\n
Optional: -p Initial probabilities based on genome wide count (g) or contig-by-contig (c) [c]
          -c Cutoff for significance [0.01]\n
Notes: maximum four numbers hard coded
       transition matrix only works for the same number consecutively and not one with repeated elements (e.g. 012023)
       https://stats.stackexchange.com/questions/26988/probability-of-finding-a-particular-sequence-of-base-pairs\n";
our($opt_c, $opt_p, $opt_q);
getopt('pqc');
die $usage unless ($opt_q);
if(!defined $opt_c) { $opt_c = 0.01; }
if(!defined $opt_p) { $opt_p = 'c'; }

# contig_gene_count = contig -> tally of genes on contig
# gene_ids = contig -> pos = quadrant
# quadrant_counts = quadrant -> gene count
my ($contig_gene_count, $gene_ids, $quadrant_counts, $gene_count) = &save_gene_count($opt_q);

# save contig -> quadrant -> gene count
# save contig -> quadrant -> consecutive count
my ($contig_quad_gene_count, $consecutive_counts) = &save_consecutive_counts_per_contig($gene_ids);


# Calculate probabilities for quadrants for whole genome
#my $q1 = ($$quadrant_counts{1} / $gene_count);
#my $q2 = ($$quadrant_counts{2} / $gene_count);
#my $q3 = ($$quadrant_counts{3} / $gene_count);
#my $q4 = ($$quadrant_counts{4} / $gene_count);
my %sig_to_line;
foreach my $contig(keys %{$consecutive_counts}) {
	foreach my $quadrant(keys %{$$consecutive_counts{$contig}}) {
		my $consec_count = $$consecutive_counts{$contig}{$quadrant};

		# ignore quadrants in a contig without any gene representatives
		next if($consec_count eq 0);

		# gene count string
		my $string;
		for(my $i=0; $i<$consec_count; $i++) {
			$string .= $quadrant;
		}

		# number of genes total
		my $contig_gene_count_total = $$contig_gene_count{$contig};
		#warn "$contig $quadrant $consec_count $string $contig_gene_count_total\n";

		# Probably not possible, but just in case string is bigger than gene count
		my $string_length = length($string);
		die "string needs to be smaller than length\n" if($string_length > $contig_gene_count_total);

		# define probability for genes in a given quadrant
		#my $q1 = ($$contig_quad_gene_count{$contig}{1} / $contig_gene_count_total);
		#my $q2 = ($$contig_quad_gene_count{$contig}{2} / $contig_gene_count_total);
		#my $q3 = ($$contig_quad_gene_count{$contig}{3} / $contig_gene_count_total);
		#my $q4 = ($$contig_quad_gene_count{$contig}{4} / $contig_gene_count_total);
		my $q;
		if($opt_p eq 'c') {
			$q = ($$contig_quad_gene_count{$contig}{$quadrant} / $contig_gene_count_total);
		} else {
			$q = ($$quadrant_counts{$quadrant} / $gene_count);
		}
		#warn "q1 = $q1, q2 = $q2, q3 = $q3 and q4 = $q4\n";
		#die;
		
		# Upper tail
		my $cumulative_probability = 0;
		CUMULATIVE_COUNT: for(my $i=0; $i<($contig_gene_count_total - $string_length); $i++) {
			my $incremental_length = ($string_length + $i);
			my $string_incremented = $string;
			for(my $h=0; $h<$i; $h++) { $string_incremented .= 'A'; }
			#print "$pretend_string\n";

			# Define transition matrix
			my @transition_matrix = &save_transition_matrix($string_incremented, $q);
			my $transition_matrix2 = Math::Matrix->new(@transition_matrix);

			# Create hash of arrays
			my %empty_matrices;
			for(my $i=0; $i<$contig_gene_count_total; $i++) {
				my @empty_array = &save_empty_matrix($string_incremented);
				my $empty_array2 = Math::Matrix->new(@empty_array);
				$empty_matrices{$i} = $empty_array2;
			}

			# let the first array be the transition_matrix
			$empty_matrices{0} = $transition_matrix2;

			# for every array (after the first), do a "matrix multiplication" (%*%) on the transition_matrix
			for(my $i=1; $i<=$contig_gene_count_total; $i++) {
				$empty_matrices{$i} = $empty_matrices{($i - 1)} -> mmul($transition_matrix2);
			}

			# Probability
			my $probability = $empty_matrices{$contig_gene_count_total}[0][$incremental_length];
			$cumulative_probability += $probability;
			#print "probability for $incremental_length length = $probability\n";

			# upper tail doesn't work (prob > 1), so don't use
			last CUMULATIVE_COUNT;
		}

		# Save
		my $save_line = "$contig\t$quadrant\t$consec_count\t$contig_gene_count_total\t$q\t$cumulative_probability\n";
		if($cumulative_probability < $opt_c) {
			$sig_to_line{$cumulative_probability} .= $save_line;
		}
		#print "$contig\t$quadrant\t$consec_count\t$contig_gene_count_total\t$q\t$cumulative_probability\n";
	}
}

# print
print "Contig\tQuadrant\tConsec_count\tTotal_count\tProb quadrant\tUpper tail probability\n";
foreach my $sig(sort { $a <=> $b } keys %sig_to_line) {
	my $lines = $sig_to_line{$sig};
	print "$lines";
}

sub save_transition_matrix {
	my ($string, $prob_letter) = @_;

	my @transition_matrix;
	my @string_parts = split //, $string;
	my $string_length = length($string);

	for(my $x=0; $x<=$string_length; $x++) {
		#warn "x = $x\n";
		for(my $y=0; $y<=$string_length; $y++) {
			#warn "y = $y\n";	
			my $x_character = $string_parts[($x - 1)];

			# default 0
			$transition_matrix[$x][$y] = 0;

			# left row
			if($y eq 0) {
				if($x < $string_length) {
					$transition_matrix[$x][$y] = 1 - $prob_letter;
				} 
			}

			# diagonal
			if(($x + 1) eq $y) {
				$transition_matrix[$x][$y] = $prob_letter;
			}

			# final position
			if(($x eq $string_length) && ($y eq $string_length)) {
				$transition_matrix[$x][$y] = 1;
			}
		
		}
	}
	return @transition_matrix;
}

sub save_empty_matrix {
	my $string = $_[0];

	my @transition_matrix;
	my $string_length = length($string);
	for(my $x=0; $x<=$string_length; $x++) {
		for(my $y=0; $y<=$string_length; $y++) {
			$transition_matrix[$x][$y] = 0;
		}
	}
	return @transition_matrix;
}

sub save_gene_count {
	my $file = $_[0];
	my (%genes, %gene_ids, %quadrant_counts);
	my $gene_count = 0;
	open my $fh, '<', $file or die "Cannot open $file : $!\n";
	while(my $line=<$fh>) {
		chomp $line;
		# ignore header
		next if($line =~ m/^scaffold\t/);

		my @bits = split /\t/, $line;
		my ($contig, $pos, $quadrant) = @bits;
		$genes{$contig}++;
		$gene_ids{$contig}{$pos} = $quadrant;
		$quadrant_counts{$quadrant}++;
		$gene_count++;
	}
	close $fh;
	return (\%genes, \%gene_ids, \%quadrant_counts, $gene_count);
}

sub save_consecutive_counts_per_contig {
	my $gene_ids = $_[0];

	# contig -> quadrant -> all consecutive count -> start = stop
	my %consecutative_quadrant_counts;
	my %contig_quad_gene_count;
	foreach my $contig(sort keys %{$gene_ids}) {
	
		my $quadrant_count = 0;
		my $previous_quadrant_found;
		my $first_position;
		my $previous_pos;

		POS: foreach my $pos(sort { $a <=> $b } keys %{$$gene_ids{$contig}}) {
			my $quadrant = $$gene_ids{$contig}{$pos};
			$contig_quad_gene_count{$contig}{$quadrant}++;

			# 1st quadrant of contig
			if(!defined $previous_quadrant_found) {
				$previous_quadrant_found = $quadrant;
				$first_position = $pos;
				$previous_pos = $pos;
				$quadrant_count++;
				next POS;
			}

			# quadrant has changed
			if($quadrant ne $previous_quadrant_found) {
			
				# save
				# contig -> quadrant -> consecutive count -> start = stop
				$consecutative_quadrant_counts{$contig}{$previous_quadrant_found}{$quadrant_count}{$first_position} = $previous_pos;

				# reset
				$previous_quadrant_found = $quadrant;
				$quadrant_count = 0;
				$first_position = $pos;
			}

			# quadrant is the same
			else {
				$quadrant_count++;
				$previous_pos = $pos;
			}
		}
	}

	# save contig -> quadrant -> largest consecutive count
	my %contig_quad_consec_count;
	foreach my $contig(sort keys %consecutative_quadrant_counts) {
		foreach my $quadrant(sort keys %{$consecutative_quadrant_counts{$contig}}) {
			QUADCOUNT: foreach my $quadrant_count(sort { $b <=> $a } keys %{$consecutative_quadrant_counts{$contig}{$quadrant}}) {
				#warn "$contig quadrant $quadrant largest = $quadrant_count\n";
				$contig_quad_consec_count{$contig}{$quadrant} = $quadrant_count;
				last QUADCOUNT;
			}
		}
	}
	return (\%contig_quad_gene_count, \%contig_quad_consec_count);
}

