#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/perl_modules";
use read_Tab;
use read_GFF;
use Data::Dumper;

### r.farrer@exeter.ac.uk

# Opening commands
my $usage = "Usage: $0 -g <GFF> > outfile\n
Parse GFF:-f\tFeature for coords (exon/CDS/gene) [gene]
          -s\tSeperator in description for gene names (\" ; etc) [;]
	  -d\tColumn that says which parent/gene it belongs to the description [0]
          -r\tRemove additional comments in column [Parent=]\n
Optional: -a\tFile1 with ids in column 0 (gene category 1)
	  -b\tFile2 with ids in column 0 (gene category 2)
	  -c\tFile3 with ids in column 0 (gene category 3)\n";
our($opt_g, $opt_f, $opt_s, $opt_d, $opt_r, $opt_a, $opt_b, $opt_c);
getopt('gfsdrabc');
die $usage unless ($opt_g);
die "$opt_g doesn't exist\n" unless (-e $opt_g);
if(!defined $opt_f) { $opt_f = 'gene'; }
if(!defined $opt_s) { $opt_s = ';'; }
if(!defined $opt_d) { $opt_d = 0; }
if(!defined $opt_r) { $opt_r = 'Parent='; }
warn "Running $0 -g $opt_g -f $opt_f -s $opt_s -d $opt_d -r $opt_r (and possibly optional files -a -b and -c)\n";

# Get additional hashes
my ($file1, $file2, $file3);
if($opt_a) { $file1 = tabfile::save_columns_to_one_hash($opt_a, 0); }
if($opt_b) { $file2 = tabfile::save_columns_to_one_hash($opt_b, 0); }
if($opt_c) { $file3 = tabfile::save_columns_to_one_hash($opt_c, 0); } 

# Save GFF and final features
my ($gff, $strand) = gfffile::gff_to_contig_parent_to_cds_hash($opt_g, $opt_f, $opt_s, $opt_d, $opt_r); # contig -> gene -> exon1 exon2
my $final_feature_per_contig = gfffile::gff_final_features($opt_g);
my $contig_to_exon_positions = &save_gff_contig_positions($gff); # contig -> start/stop = gene 

# Calculate 5' and 3' distance to next gene
my $count = 0;
my $count_non_terminal_entries = 0;
print "contig\tfeature\tdesc\tstrand\tupstream_distance\tdownstream_distance\tstart_pos\tuplog10\tdownlog10\n";
CONTIG: foreach my $contig(sort keys %{$gff}) {
	GENES: foreach my $gene(sort keys %{$$gff{$contig}}) {

		# 1st gene, to check on ID
		if($count eq 0) {
			warn "1st gene = $gene\n";
			$count = 1;
		}

		my $start_stop = $$gff{$contig}{$gene};
		my @start_and_stops = split / /, $start_stop;

		#warn "looking at $contig $gene with @start_and_stops...\n";
		my $start_pos_saved;
		EXONS: foreach my $start_and_stop(@start_and_stops) {
			my @exon_bounderies = split /-/, $start_and_stop;
			my $start_pos = $exon_bounderies[0];
			my $stop_pos = $exon_bounderies[1];
			my $strand = $$strand{$gene};

			# save start position (new)
			if(!defined $start_pos_saved) { $start_pos_saved = $start_pos; }

			# description
			my $desc;
			if(defined($$file1{$gene})) { $desc .= "cat1,"; }
			if(defined($$file2{$gene})) { $desc .= "cat2,"; }
			if(defined($$file3{$gene})) { $desc .= "cat3,"; }
			if(!defined $desc) { $desc = 'NA'; }

			# Entry to the right
			# The ATG of the gene is from 3'->5'. 
			# If gene is + and to the right = upstream / 5'
			# if gene is - and to the right = downstream / 3'
			# if gene is + and to the left = downstream / 3'
			# if gene is - and to the left = upstream / 5'
			my ($upstream_distance, $downstream_distance) = (0, 0);
			
			# Right
			THREEPRIME: for(my $i=$stop_pos; $i<$$final_feature_per_contig{$contig}; $i++) {

				# Found something that's not itself
				if((defined $$contig_to_exon_positions{$contig}{$i}) && ($$contig_to_exon_positions{$contig}{$i} ne $gene)) {	
					my $distance_to_right = ($i - $stop_pos);
					if($strand eq '+') { $upstream_distance = $distance_to_right; }
					else { $downstream_distance = $distance_to_right; }
					last THREEPRIME;
				}
			}

			# Left
			FIVEPRIME: for(my $i=($start_pos - 1); $i > 1; $i--) {

				# Found something that's not itself
				if((defined $$contig_to_exon_positions{$contig}{$i}) && ($$contig_to_exon_positions{$contig}{$i} ne $gene)) {	
					my $distance_to_left = ($start_pos - $i);
					if($strand eq '+') { $downstream_distance = $distance_to_left; }
					else { $upstream_distance = $distance_to_left; }
					last FIVEPRIME;
				}
			}

			# no other genes on this contig or flanking
			#warn "upstream = $upstream_distance and downstream = $downstream_distance\n";
			next GENES if(($upstream_distance eq 0) || ($downstream_distance eq 0));
			$count_non_terminal_entries++;

			# Log
			my ($log_up, $log_down) = (0, 0);
			if($upstream_distance ne 0) { $log_up = &log10($upstream_distance); }
			if($downstream_distance ne 0) { $log_down = &log10($downstream_distance); }

			# Print
			print "$contig\t$gene\t$desc\t$strand\t$upstream_distance\t$downstream_distance\t$start_pos_saved\t$log_up\t$log_down\n";
		}
	}
}
warn "$count_non_terminal_entries features that are not terminal entries\n"; 

sub log10 {
	my $n = shift;
	return (log($n) / log(10));
}

sub save_gff_contig_positions {
	my $gff = $_[0];
	my %exon_positions;
	warn "save_gff_contig_positions...\n";
	CONTIG: foreach my $contig(sort keys %{$gff}) {
		GENES: foreach my $gene(keys %{$$gff{$contig}}) {
			my $start_stop = $$gff{$contig}{$gene};
			my @start_and_stops = split / /, $start_stop;
			#warn "looking at $contig $gene with @start_and_stops...\n";
			EXONS: foreach my $start_and_stop(@start_and_stops) {
				my @exon_bounderies = split /-/, $start_and_stop;
				my $start_pos = $exon_bounderies[0];
				my $stop_pos = $exon_bounderies[1];

				# Save
				#warn "Saving $contig $start_pos and $stop_pos.\n";
				$exon_positions{$contig}{$start_pos} = $gene;
				$exon_positions{$contig}{$stop_pos} = $gene;
			}
		}
	}
	return \%exon_positions;
}
