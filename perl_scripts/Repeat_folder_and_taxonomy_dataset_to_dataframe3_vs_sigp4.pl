#!/usr/bin/perl -w
use strict;
use Getopt::Std;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -r <folder of assemblies that have been repeatmodelled, repeatmasked and a summary made + .SigP4> -t <taxanomy_data_for_accessions.tab> -g <genome_lengths.tab>\n
Notes: Makes 2 outfile dataframes\n";
our($opt_r, $opt_t, $opt_g);
getopt('rtg');
die $usage unless ($opt_r && $opt_t && $opt_g);
foreach($opt_r, $opt_t, $opt_g) { die "Cannot open $_ : $!" unless (-e $_); }

# Output files
my $outfile1 = "./RF_12345_Repeat_vs_secreted_flanking.tab";
my $outfile2 = "./RF_12345_Repeat_vs_secreted_not_flanking.tab";
my $outfile3 = "./RF_12345_Repeat_vs_secreted_flanking_summary.tab";
my $outfile4 = "./RF_12345_Repeat_vs_secreted_not_flanking_summary.tab";

my $outfile1_upstream = "./RF_12345_Repeat_vs_secreted_flanking_upstream.tab";
my $outfile2_upstream = "./RF_12345_Repeat_vs_secreted_not_flanking_upstream.tab";
my $outfile3_upstream = "./RF_12345_Repeat_vs_secreted_flanking_summary_upstream.tab";
my $outfile4_upstream = "./RF_12345_Repeat_vs_secreted_not_flanking_summary_upstream.tab";

my $outfile1_downstream = "./RF_12345_Repeat_vs_secreted_flanking_downstream.tab";
my $outfile2_downstream = "./RF_12345_Repeat_vs_secreted_not_flanking_downstream.tab";
my $outfile3_downstream = "./RF_12345_Repeat_vs_secreted_flanking_downstream_summary.tab";
my $outfile4_downstream = "./RF_12345_Repeat_vs_secreted_not_flanking_downstream_summary.tab";

my $outfile_HgT = "./RF_12345_Repeat_vs_secreted_HgT_info.tab";
my $outfile_HgT2 = "./RF_12345_Repeat_vs_secreted_HgT_info_summary.tab";

my $outfile1_upstream_distances = "./RF_12345_Repeat_vs_secreted_flanking_upstream_distances.tab";
my $outfile1_downstream_distances = "./RF_12345_Repeat_vs_secreted_flanking_downstream_distances.tab";
my $outfile2_upstream_distances = "./RF_12345_Repeat_vs_secreted_not_flanking_upstream_distances.tab";
my $outfile2_downstream_distances = "./RF_12345_Repeat_vs_secreted_not_flanking_downstream_distances.tab";

#my $outfile_number_of_repeats_excluded_from_analysis = "./RF_12345_Repeat_vs_secreted_repeats_excluded.tab";

my $outfile_repeats_flanking_secreted = "./RF_12345_Repeat_vs_secreted_repeats_flanking_secreted_tally.tab";
my $outfile_repeats_not_flanking_secreted = "./RF_12345_Repeat_vs_secreted_repeats_not_flanking_secreted_tally.tab";

# Check what has already computed
my $assemblies_computed = &check_which_dataframes_have_been_computed($outfile1);

# All data
# assembly -> 'kingdom' -> name
# assembly -> species_name -> name
my $assembly_to_taxonomy = &save_taxonomy_data($opt_t);

# Assembly -> length = length
# Assembly -> n50 = n50
my $assembly_metrics = &save_genome_lengths($opt_g);

# Save all repeat folders
warn "Save each assembly with repeats in $opt_r...\n";
my @assemblies_with_repeats_info = <$opt_r/*>;
ASSEMBLY: foreach my $assembly_dir(@assemblies_with_repeats_info) {

	# split e.g. 9.Repeats_all/finished//GCA_000001985.1
	my @assembly_dir_parts = split /\//, $assembly_dir;
	my $assembly = $assembly_dir_parts[-1];
	#warn "assembly = $assembly\n";
	die "No assembly found for $assembly (from folder) in taxonomy file\n" if(!defined $$assembly_to_taxonomy{$assembly});
	die "no kingdom for $assembly (from folder) in taxonomy file\n" if(!defined $$assembly_to_taxonomy{$assembly}{'kingdom'});
	die "no species name for $assembly (from folder) in taxonomy file\n" if(!defined $$assembly_to_taxonomy{$assembly}{'species_name'});
	my $kingdom = $$assembly_to_taxonomy{$assembly}{'kingdom'};
	my $species = $$assembly_to_taxonomy{$assembly}{'species_name'};

	# Check if the files already exist. if so, lets only re-calculate those not evaluated, as it's extremely slow. Then append new ones on the old file
	if(defined $$assemblies_computed{$assembly}) {
		warn "Skipping: $assembly - already computed based on presence in $outfile1\n";
		next ASSEMBLY;
	}

	# assembly metrics
	die "no genome length for $assembly\n" if(!defined $$assembly_metrics{$assembly}{'length'});
	my $genome_length = $$assembly_metrics{$assembly}{'length'};

	# Assembly -> gene = SignalP4 (1 or no)
	my @secreted_info = <$assembly_dir/*.SigP.clean>;
	die "no repeat info found in $assembly_dir\n" if(!defined $secreted_info[0]);
	die "more than 1 SigP4 found in $assembly_dir\n" if(defined $secreted_info[1]);
	my $sigp4_file = $secreted_info[0];
	my $genes_with_sigp4 = &save_gene_to_sigp4($sigp4_file, $assembly);

	# save repeat locations (assembly -> contig -> "$query_start-$query_stop" = repeat class)
	my @repeat_info = <$assembly_dir/*.out>;
	die "no repeat info found in $assembly_dir\n" if(!defined $repeat_info[0]);
	die "more than 1 repeat info found in $assembly_dir\n" if(defined $repeat_info[1]);
	my $repeat_file = $repeat_info[0];
	my ($repeat_summary_info, $repeat_summary_names) = &save_repeat_info($repeat_file, $assembly);

	# save coords to make sense of sigp results (Assembly -> contig -> gene start or stop = gene_id
	# coords_final_feature = Assembly -> contig = final base with annotation
	# strand =  gene = strand
	my @coords_info = <$assembly_dir/*.coords>;
	die "no repeat info found in $assembly_dir\n" if(!defined $coords_info[0]);
	die "more than 1 .coords found in $assembly_dir\n" if(defined $coords_info[1]);
	my $coords_file = $coords_info[0];
	my ($coords_summary_info, $coords_final_feature, $coords_strand) = &save_coords_info($coords_file, $assembly);

	# dataframe info
	my %dataframe_info;
	$dataframe_info{$assembly}{'kingdom'} = $kingdom;
	$dataframe_info{$assembly}{'species'} = $species;
	$dataframe_info{$assembly}{'genome_length'} = $genome_length;

	# Track which genes were already counted as flanking a repeat.
	# Each gene should be only counted once, even if it flanks multiple repeats. Just take the 1st repeat whatever it is!
	my %genes_already_assigned_to_repeats;

	# flanking a secreted? (either side)
	my %repeat_type_info_secreted; # assembly -> repeat_name -> nt and percent
	my %repeat_type_info_not_secreted; # assembly -> repeat_name -> nt and percent
	my %repeat_type_info_secreted_summary; # same but condensed (e.g., all DNA transposons joined)
	my %repeat_type_info_not_secreted_summary; # same but condensed (e.g., all DNA transposons joined)

	# flanking a secreted upstream (maybe on 1 side or both sides of the repeat)
	my %repeat_type_info_secreted_upstream; # assembly -> repeat_name -> nt and percent
	my %repeat_type_info_not_secreted_upstream; # assembly -> repeat_name -> nt and percent
	my %repeat_type_info_secreted_summary_upstream; # same but condensed (e.g., all DNA transposons joined)
	my %repeat_type_info_not_secreted_summary_upstream; # same but condensed (e.g., all DNA transposons joined)

	# flanking a secreted downstream (maybe on 1 side or both sides of the repeat)
	my %repeat_type_info_secreted_downstream; # assembly -> repeat_name -> nt and percent
	my %repeat_type_info_not_secreted_downstream; # assembly -> repeat_name -> nt and percent
	my %repeat_type_info_secreted_summary_downstream; # same but condensed (e.g., all DNA transposons joined)
	my %repeat_type_info_not_secreted_summary_downstream; # same but condensed (e.g., all DNA transposons joined)

	# Hypergeometric test lines (Assembly -> repeat = line)
	my @HgT_info;
	my %HgT_info2; # q1.repeat_flanking_secreted_gene, q2.repeat_flanking_upstream_secreted_gene, q3.repeat_flanking_downstream_secreted_gene
	               # m.total_secreted_genes
	               # n.total_non_secreted_genes
	               # k.total_number_of_flanking_genes

	my %upstream_distances_flanking_secreted;
	my %downstream_distances_flanking_secreted;
	my %upstream_distances_not_flanking_secreted;
	my %downstream_distances_not_flanking_secreted;

	#my %ignored_repeat_classes;
	#my $unaccounted_repeats = 0;

	my %repeat_flanking_secreted_tallies;
	my %repeat_not_flanking_secreted_tallies;

	# now to determine flanking repeats present between genes?
	warn "Determining which genes are flanking repeats...\n";
	foreach my $contig(sort keys %{$$repeat_summary_info{$assembly}}) {
		#warn "\tcontig $contig (unaccounted repeats = $unaccounted_repeats)\n";
		warn "\tcontig $contig\n";


		REPEATS: foreach my $repeat_bounderies(sort keys %{$$repeat_summary_info{$assembly}{$contig}}) {
			my $repeat_name = $$repeat_summary_names{$assembly}{$contig}{$repeat_bounderies};
			my $repeat_class = $$repeat_summary_info{$assembly}{$contig}{$repeat_bounderies};
			my @repeat_class_split = split /\//, $repeat_class;
			my $repeat_class_type = $repeat_class_split[0]; 
			my @repeat_bounderies_split = split /-/, $repeat_bounderies;
			my ($repeat_start, $repeat_end) = @repeat_bounderies_split;
			my $repeat_length = ($repeat_end - $repeat_start);
			#warn "repeat $repeat_class bounderies = $contig $repeat_start $repeat_end\n";

			# The ATG of the gene is from 3'->5'. 
			# If gene is + and to the right = upstream / 5'
			# if gene is - and to the right = downstream / 3'
			# if gene is + and to the left = downstream / 3'
			# if gene is - and to the left = upstream / 5'
			my ($upstream_distance, $downstream_distance) = (0, 0);

			# go right till i find a gene, and see if that gene is secreted
			my $found_gene_right = 0;
			my $found_gene_right_id = 'None';
			my $found_gene_right_strand = '?';
			my $found_gene_right_distance = 0;
			my $found_gene_right_secreted = 'N';
			if(defined $$coords_final_feature{$assembly}{$contig}) {
				RIGHT: for(my $i=$repeat_end; $i < $$coords_final_feature{$assembly}{$contig}; $i++) {
					if(defined $$coords_summary_info{$assembly}{$contig}{$i}) {
						$found_gene_right = 1;
						$found_gene_right_id = $$coords_summary_info{$assembly}{$contig}{$i};
						$found_gene_right_strand = $$coords_strand{$found_gene_right_id};
						$found_gene_right_distance = ($i - $repeat_end);
						if(defined $$genes_with_sigp4{$assembly}{$found_gene_right_id}) { $found_gene_right_secreted = 'Y'; }
						#warn "FOUND the gene flanking to the right: $i = $found_gene_right_id\n";
						last RIGHT;
					}
				}
			}

			# go left till i find a gene, and see if that gene is secreted
			my $found_gene_left = 0;
			my $found_gene_left_id = 'None';
			my $found_gene_left_strand = '?';
			my $found_gene_left_distance = 0;
			my $found_gene_left_secreted = 'N';
			LEFT: for(my $i=$repeat_start; $i > 0; $i--) {
				if(defined $$coords_summary_info{$assembly}{$contig}{$i}) {
					$found_gene_left = 1;
					$found_gene_left_id = $$coords_summary_info{$assembly}{$contig}{$i};
					$found_gene_left_strand = $$coords_strand{$found_gene_left_id};
					$found_gene_left_distance = ($repeat_start - $i);
					if(defined $$genes_with_sigp4{$assembly}{$found_gene_left_id}) { $found_gene_left_secreted = 'Y'; }
					#warn "FOUND the gene flanking to the left: $i = $found_gene_left_id\n";
					last LEFT;
				}
			}

			# flanking a gene encoding a secreted protein
			if(($found_gene_right_secreted eq 'Y') || ($found_gene_left_secreted eq 'Y')) {
				$repeat_type_info_secreted{$repeat_class}{$assembly} += $repeat_length;
				$repeat_type_info_secreted_summary{$repeat_class_type}{$assembly} += $repeat_length;
				$repeat_flanking_secreted_tallies{$repeat_class}++;
			} else {
				$repeat_type_info_not_secreted{$repeat_class}{$assembly} += $repeat_length;
				$repeat_type_info_not_secreted_summary{$repeat_class_type}{$assembly} += $repeat_length;
				$repeat_not_flanking_secreted_tallies{$repeat_class}++;
			}

			# flanking a secreted upstream (maybe on 1 side or both sides of the repeat)
			if((($found_gene_right_secreted eq 'Y') && ($found_gene_right_strand eq '+')) || 
			   (($found_gene_left_secreted eq 'Y') && ($found_gene_left_strand eq '-'))) {
				$repeat_type_info_secreted_upstream{$repeat_class}{$assembly} += $repeat_length;
				$repeat_type_info_secreted_summary_upstream{$repeat_class_type}{$assembly} += $repeat_length;
			} else {
				$repeat_type_info_not_secreted_upstream{$repeat_class}{$assembly} += $repeat_length;
				$repeat_type_info_not_secreted_summary_upstream{$repeat_class_type}{$assembly} += $repeat_length;
			}

			# flanking a secreted downstream (maybe on 1 side or both sides of the repeat)
			if((($found_gene_right_secreted eq 'Y') && ($found_gene_right_strand eq '-')) || 
			   (($found_gene_left_secreted eq 'Y') && ($found_gene_left_strand eq '+'))) {
				$repeat_type_info_secreted_downstream{$repeat_class}{$assembly} += $repeat_length;
				$repeat_type_info_secreted_summary_downstream{$repeat_class_type}{$assembly} += $repeat_length;
			} else {
				$repeat_type_info_not_secreted_downstream{$repeat_class}{$assembly} += $repeat_length;
				$repeat_type_info_not_secreted_summary_downstream{$repeat_class_type}{$assembly} += $repeat_length;
			}

			# Distances
			if(($found_gene_right_secreted eq 'Y') && ($found_gene_right_strand eq '+')) {
				$upstream_distances_flanking_secreted{$found_gene_right_distance}++;
			} else {
				$upstream_distances_not_flanking_secreted{$found_gene_right_distance}++;
			}
			if(($found_gene_left_secreted eq 'Y') && ($found_gene_left_strand eq '-')) {
				$upstream_distances_flanking_secreted{$found_gene_left_distance}++;
			} else {
				$upstream_distances_not_flanking_secreted{$found_gene_left_distance}++;
			}

			if(($found_gene_right_secreted eq 'Y') && ($found_gene_right_strand eq '-')) {
				$downstream_distances_flanking_secreted{$found_gene_right_distance}++;
			} else {
				$downstream_distances_not_flanking_secreted{$found_gene_right_distance}++;
			}
			if(($found_gene_left_secreted eq 'Y') && ($found_gene_left_strand eq '+')) {
				$downstream_distances_flanking_secreted{$found_gene_left_distance}++;
			} else {
				$downstream_distances_not_flanking_secreted{$found_gene_left_distance}++;
			}

			### DATAFRAME FOR HYPERGEOMETRIC TESTS

			### each gene set only counted once (mention excluding these for HgT only)
			if((defined $genes_already_assigned_to_repeats{$found_gene_right_id}{$found_gene_left_id}) ||
				(defined $genes_already_assigned_to_repeats{$found_gene_left_id}{$found_gene_right_id})) {
					#$unaccounted_repeats++;
					#$ignored_repeat_classes{$repeat_class_type}++;
					next REPEATS;
			}

			# Assembly | contig | repeat_id | left_gene | right_gene | left_strand | right_strand | left_SP | right_SP | category
			#1) 2 genes upstream of it (either side, genes facing away)
			#2) 2 genes downstream of it (either side, genes facing towards)
			#3) 1 gene upstream and 1 gene downstream
			#4) 1 gene upstream only, as repeat is at the end of the contig (no further annotation)
			#5) 1 gene downstream only, same as above
			#6) nothing, as repeat is on contig without any genes

			my $category;
			if(($found_gene_right_strand eq '+') && ($found_gene_left_strand eq '-')) { 
				$category = '1.both_genes_upstream'; 
				$HgT_info2{'k.total_number_of_flanking_genes'} += 2;
				if($found_gene_right_secreted eq 'Y') { $HgT_info2{'q2.repeat_flanking_upstream_secreted_gene'}++; }
				if($found_gene_left_secreted eq 'Y') { $HgT_info2{'q2.repeat_flanking_upstream_secreted_gene'}++; }
			}
			elsif(($found_gene_right_strand eq '-') && ($found_gene_left_strand eq '+')) { 
				$category = '2.both_genes_downstream'; 
				$HgT_info2{'k.total_number_of_flanking_genes'} += 2;
				if($found_gene_right_secreted eq 'Y') { $HgT_info2{'q3.repeat_flanking_downstream_secreted_gene'}++; }
				if($found_gene_left_secreted eq 'Y') { $HgT_info2{'q3.repeat_flanking_downstream_secreted_gene'}++; }
			}
			elsif(($found_gene_right_strand eq '+') && ($found_gene_left_strand eq '+')) { 
				$category = '3.one_upstream_one_downstream'; 
				$HgT_info2{'k.total_number_of_flanking_genes'} += 2;
				if($found_gene_right_secreted eq 'Y') { $HgT_info2{'q2.repeat_flanking_upstream_secreted_gene'}++; }
				if($found_gene_left_secreted eq 'Y') { $HgT_info2{'q3.repeat_flanking_downstream_secreted_gene'}++; }
			}
			elsif(($found_gene_right_strand eq '-') && ($found_gene_left_strand eq '-')) { 
				$category = '3.one_upstream_one_downstream'; 
				$HgT_info2{'k.total_number_of_flanking_genes'} += 2;
				if($found_gene_right_secreted eq 'Y') { $HgT_info2{'q3.repeat_flanking_downstream_secreted_gene'}++; }
				if($found_gene_left_secreted eq 'Y') { $HgT_info2{'q2.repeat_flanking_upstream_secreted_gene'}++; }
			}
			elsif(($found_gene_right_strand eq '+') && ($found_gene_left_id eq 'None')) { 
				$category = '4.one_upstream'; 
				$HgT_info2{'k.total_number_of_flanking_genes'}++;
				if($found_gene_right_secreted eq 'Y') { $HgT_info2{'q2.repeat_flanking_upstream_secreted_gene'}++; }
			}
			elsif(($found_gene_left_strand eq '-') && ($found_gene_right_id eq 'None')) { 
				$category = '4.one_upstream'; 
				$HgT_info2{'k.total_number_of_flanking_genes'}++;
				if($found_gene_left_secreted eq 'Y') { $HgT_info2{'q2.repeat_flanking_upstream_secreted_gene'}++; }
			}
			elsif(($found_gene_right_strand eq '-') && ($found_gene_left_id eq 'None')) { 
				$category = '5.one_downstream'; 
				$HgT_info2{'k.total_number_of_flanking_genes'}++;
				if($found_gene_right_secreted eq 'Y') { $HgT_info2{'q3.repeat_flanking_downstream_secreted_gene'}++; }
			}
			elsif(($found_gene_left_strand eq '+') && ($found_gene_right_id eq 'None')) { 
				$category = '5.one_downstream'; 
				$HgT_info2{'k.total_number_of_flanking_genes'}++;
				if($found_gene_left_secreted eq 'Y') { $HgT_info2{'q3.repeat_flanking_downstream_secreted_gene'}++; }
			}
			elsif(($found_gene_left_id eq 'None') && ($found_gene_right_id eq 'None')) { $category = '6.no_genes_found'; }
			else {
				warn "WARNING: No category found for left ($found_gene_left_id, $found_gene_left_strand) and right ($found_gene_right_id, $found_gene_right_strand)\n";
				$category = '7.unknown';
			}

			# other HgT info2
			# q1.repeat_flanking_secreted_gene, q2.repeat_flanking_upstream_secreted_gene, q3.repeat_flanking_downstream_secreted_gene
	        # m.total_secreted_genes
	        # n.total_non_secreted_genes
	        # k.total_number_of_flanking_genes
			if($found_gene_right_secreted eq 'Y') { 
				$HgT_info2{'q1.repeat_flanking_secreted_gene'}++;
				$HgT_info2{'m.total_secreted_genes'}++; 
			}
			else { $HgT_info2{'n.total_non_secreted_genes'}++; }

			if($found_gene_left_secreted eq 'Y') { 
				$HgT_info2{'q1.repeat_flanking_secreted_gene'}++;
				$HgT_info2{'m.total_secreted_genes'}++; 
			}
			else { $HgT_info2{'n.total_non_secreted_genes'}++; }

			

			# Given we're ignoring so many repeats (millions!) - the repeat class and name isn't very meaningful here. So, consider it representing all those not included.
			my $HgT_line = "$assembly\t$contig\t$repeat_class\t$repeat_name\t$found_gene_left_id\t$found_gene_right_id\t$found_gene_left_strand\t$found_gene_right_strand\t$found_gene_left_secreted\t$found_gene_right_secreted\t$category\n";
			push @HgT_info, $HgT_line;


			# Record genes seen (Avoid multiple repeats in a given intergenic space being counted twice)
			$genes_already_assigned_to_repeats{$found_gene_right_id}{$found_gene_left_id} = 1;
			$genes_already_assigned_to_repeats{$found_gene_left_id}{$found_gene_right_id} = 1;
			#die "end here \n";
		}
	}
	#die "end here";

	### flanking a secreted? (either side)
	# print headers (if not exists yet)
	&print_dataframes_header(\%repeat_type_info_secreted, $outfile1);
	&print_dataframes_header(\%repeat_type_info_not_secreted, $outfile2);
	&print_dataframes_header(\%repeat_type_info_secreted_summary, $outfile3);
	&print_dataframes_header(\%repeat_type_info_not_secreted_summary, $outfile4);

	# print dataframes (for plots)
	&print_dataframes(\%dataframe_info, \%repeat_type_info_secreted, $outfile1);
	&print_dataframes(\%dataframe_info, \%repeat_type_info_not_secreted, $outfile2);
	&print_dataframes(\%dataframe_info, \%repeat_type_info_secreted_summary, $outfile3);
	&print_dataframes(\%dataframe_info, \%repeat_type_info_not_secreted_summary, $outfile4);

	### flanking a secreted upstream (maybe on 1 side or both sides of the repeat)
	# print headers (if not exists yet)
	&print_dataframes_header(\%repeat_type_info_secreted_upstream, $outfile1_upstream);
	&print_dataframes_header(\%repeat_type_info_not_secreted_upstream, $outfile2_upstream);
	&print_dataframes_header(\%repeat_type_info_secreted_summary_upstream, $outfile3_upstream);
	&print_dataframes_header(\%repeat_type_info_not_secreted_summary_upstream, $outfile4_upstream);

	# print dataframes (for plots)
	&print_dataframes(\%dataframe_info, \%repeat_type_info_secreted_upstream, $outfile1_upstream);
	&print_dataframes(\%dataframe_info, \%repeat_type_info_not_secreted_upstream, $outfile2_upstream);
	&print_dataframes(\%dataframe_info, \%repeat_type_info_secreted_summary_upstream, $outfile3_upstream);
	&print_dataframes(\%dataframe_info, \%repeat_type_info_not_secreted_summary_upstream, $outfile4_upstream);

	### flanking a secreted downstream (maybe on 1 side or both sides of the repeat)
	# print headers (if not exists yet)
	&print_dataframes_header(\%repeat_type_info_secreted_downstream, $outfile1_downstream);
	&print_dataframes_header(\%repeat_type_info_not_secreted_downstream, $outfile2_downstream);
	&print_dataframes_header(\%repeat_type_info_secreted_summary_downstream, $outfile3_downstream);
	&print_dataframes_header(\%repeat_type_info_not_secreted_summary_downstream, $outfile4_downstream);

	# print dataframes (for plots)
	&print_dataframes(\%dataframe_info, \%repeat_type_info_secreted_downstream, $outfile1_downstream);
	&print_dataframes(\%dataframe_info, \%repeat_type_info_not_secreted_downstream, $outfile2_downstream);
	&print_dataframes(\%dataframe_info, \%repeat_type_info_secreted_summary_downstream, $outfile3_downstream);
	&print_dataframes(\%dataframe_info, \%repeat_type_info_not_secreted_summary_downstream, $outfile4_downstream);

	# print HgT info
	&print_hgt_dataframes_header($outfile_HgT);
	open my $ofh, '>>', $outfile_HgT or die "error: Cannot open $outfile_HgT : $!";
	foreach my $HgT_line(@HgT_info) {
		print $ofh "$HgT_line";
	}
	close $ofh;

	# print HgT info2
	if(!defined $HgT_info2{'q1.repeat_flanking_secreted_gene'}) { $HgT_info2{'q1.repeat_flanking_secreted_gene'} = 0; }
	if(!defined $HgT_info2{'q2.repeat_flanking_upstream_secreted_gene'}) { $HgT_info2{'q2.repeat_flanking_upstream_secreted_gene'} = 0; }
	if(!defined $HgT_info2{'q3.repeat_flanking_downstream_secreted_gene'}) { $HgT_info2{'q3.repeat_flanking_downstream_secreted_gene'} = 0; }
	if(!defined $HgT_info2{'m.total_secreted_genes'}) { $HgT_info2{'m.total_secreted_genes'} = 0; }
	if(!defined $HgT_info2{'k.total_number_of_flanking_genes'}) { $HgT_info2{'k.total_number_of_flanking_genes'} = 0; }
	open my $ofh2, '>>', $outfile_HgT2 or die "error: Cannot open $outfile_HgT2 : $!";
	print $ofh2 "### $assembly\n";
	#print $ofh2 "# flanking secreted\n"; # DOESNT MAKE SENSE BECAUSE EVERY SECRETED GENE IS FLANKED BY A REPEAT!
	#print $ofh2 "phyper($HgT_info2{'q1.repeat_flanking_secreted_gene'}, $HgT_info2{'m.total_secreted_genes'}, $HgT_info2{'n.total_non_secreted_genes'}, $HgT_info2{'k.total_number_of_flanking_genes'}, lower.tail = FALSE)\n";
	print $ofh2 "# flanking secreted (upstream)\n";
	print $ofh2 "phyper($HgT_info2{'q2.repeat_flanking_upstream_secreted_gene'}, $HgT_info2{'m.total_secreted_genes'}, $HgT_info2{'n.total_non_secreted_genes'}, $HgT_info2{'k.total_number_of_flanking_genes'}, lower.tail = FALSE)\n";
	print $ofh2 "# flanking secreted (downstream)\n";
	print $ofh2 "phyper($HgT_info2{'q3.repeat_flanking_downstream_secreted_gene'}, $HgT_info2{'m.total_secreted_genes'}, $HgT_info2{'n.total_non_secreted_genes'}, $HgT_info2{'k.total_number_of_flanking_genes'}, lower.tail = FALSE)\n";
	print $ofh2 "\n";
	close $ofh2;

	# print distances
	&print_distances_dataframes_header($outfile1_upstream_distances);
	&print_distances_dataframes_header($outfile1_downstream_distances);
	&print_distances_dataframes_header($outfile2_upstream_distances);
	&print_distances_dataframes_header($outfile2_downstream_distances);

	&print_distances_dataframe($outfile1_upstream_distances, \%upstream_distances_flanking_secreted, $assembly);
	&print_distances_dataframe($outfile1_downstream_distances, \%downstream_distances_flanking_secreted, $assembly);
	&print_distances_dataframe($outfile2_upstream_distances, \%upstream_distances_not_flanking_secreted, $assembly);
	&print_distances_dataframe($outfile2_downstream_distances, \%downstream_distances_not_flanking_secreted, $assembly);

	# ignored repeats
	#open my $ofh3, '>>', $outfile_number_of_repeats_excluded_from_analysis or die "error: Cannot open $outfile_number_of_repeats_excluded_from_analysis : $!";
	#foreach my $class(sort keys %ignored_repeat_classes) {
	#	my $tally = $ignored_repeat_classes{$class};
	#	print $ofh3 "$assembly\t$class\t$tally\n";
	#}
	#close $ofh3;

	# repeat tallies (numbers of them, rather than nt covered)
	&print_tallies_dataframes_header($outfile_repeats_flanking_secreted);
	&print_tallies_dataframes_header($outfile_repeats_not_flanking_secreted);
	&print_distances_dataframe($outfile_repeats_flanking_secreted, \%repeat_flanking_secreted_tallies, $assembly);
	&print_distances_dataframe($outfile_repeats_not_flanking_secreted, \%repeat_not_flanking_secreted_tallies, $assembly);

}

sub print_tallies_dataframes_header {
	my ($outfile) = @_;

	# Don't print header if file already exists
	return if(-e $outfile);

	# open out file handles
	open my $ofh, '>', $outfile or die "error: Cannot open $outfile : $!";

	# header
	print $ofh "Assembly\repeat_family\ttally\n";
	close $ofh;
	return;
}

sub print_distances_dataframes_header {
	my ($outfile) = @_;

	# Don't print header if file already exists
	return if(-e $outfile);

	# open out file handles
	open my $ofh, '>', $outfile or die "error: Cannot open $outfile : $!";

	# header
	print $ofh "Assembly\tdistance\ttally\n";
	close $ofh;
	return;
}

sub print_distances_dataframe {
	my ($outfile, $distances, $assembly) = @_;

	# open out file handles (Append)
	open my $ofh, '>>', $outfile or die "error: Cannot open $outfile : $!";
	foreach my $distance(sort keys %{$distances}) {
		my $tally = $$distances{$distance};
		print $ofh "$assembly\t$distance\t$tally\n";
	}
	close $ofh;
	return;
}

sub print_hgt_dataframes_header {
	my ($outfile) = @_;

	# Don't print header if file already exists
	return if(-e $outfile);

	# open out file handles
	open my $ofh, '>', $outfile or die "error: Cannot open $outfile : $!";

	# header
	print $ofh "Assembly\tcontig\tREPRESENTING_repeat_class\tREPRESENTING_repeat_id\tleft_gene\tright_gene\tleft_strand\tright_strand\tleft_SP\tright_SP\tcategory\n";
	close $ofh;
	return;
}

sub print_dataframes {
	my ($dataframe_info, $repeat_df, $outfile) = @_;

	# open out file handles (Append)
	open my $ofh, '>>', $outfile or die "error: Cannot open $outfile : $!";

	# data
	foreach my $assembly(sort keys %{$dataframe_info}) {
		my $kingdom = $$dataframe_info{$assembly}{'kingdom'};
		my $species = $$dataframe_info{$assembly}{'species'};
		my $genome_length = $$dataframe_info{$assembly}{'genome_length'};
		print $ofh "$assembly\t$kingdom\t$species\t$genome_length";
		foreach my $repeat(sort keys %{$repeat_df}) {
			my $tally = 0;
			if(defined $$repeat_df{$repeat}{$assembly}) { $tally = $$repeat_df{$repeat}{$assembly}; }
			print $ofh "\t$tally";
		}
		print $ofh "\n";
	}
	close $ofh;
	return;
}

sub print_dataframes_header {
	my ($repeat_df, $outfile) = @_;

	# Don't print header if file already exists (DO BECAUSE IT MIGHT DIFFER BETWEEN DATASETS)
	#return if(-e $outfile);

	# open out file handles (append)
	open my $ofh, '>>', $outfile or die "error: Cannot open $outfile : $!";

	# header
	print $ofh "Assembly\tKingdom\tSpecies\tGenome_length";
	foreach my $repeat(sort keys %{$repeat_df}) {
		print $ofh "\t$repeat";
	}
	print $ofh "\n";
	close $ofh;
	return;
}

sub check_which_dataframes_have_been_computed {
	my ($infile) = @_;

	my %assemblies;

	# nothing yet!
	return if(! -e $infile);

	# open in file
	open my $fh, '<', $infile or die "error: Cannot open $infile : $!";
	while(my $line=<$fh>) {
		chomp $line;

		# header
		next if($line =~ m/^Assembly\tKingdom\tSpecies/);

		my @bits = split /\t/, $line;
		$assemblies{$bits[0]} = 1;
	}

	return \%assemblies;
}

sub save_taxonomy_data {
	my $file = $_[0];
	my %data;

	# save taxonomy data E.g.
	# Assembly_accession	TaxID	Species_name	Strain	clade	class	genus	kingdom
	# GCA_000001215.4	7227	Drosophila melanogaster		Opisthokonta	Insecta	Drosophila	Metazoa
	# GCA_000001735.2	3702	Arabidopsis thaliana	ecotype=Columbia	Embryophyta	Magnoliopsida	Arabidopsis	Viridiplantae
	open my $fh, '<', $file or die "Cannot open $file : $!\n";
	while(my $line=<$fh>) {
		chomp $line;

		# ignore header
		next if($line =~ m/^Assembly_accession/);

		my @bits = split /\t/, $line;
		my ($assembly, $taxID, $species_name, $strain, $clade, $class, $genus, $kingdom) = @bits;
		$data{$assembly}{'kingdom'} = $kingdom;
		$data{$assembly}{'species_name'} = $species_name;
	}
	return \%data;
}

sub save_genome_lengths {
	my $file = $_[0];
	my %data;

	# save genome lengths and N50 E.g.
	# Eukaryota/GCA_000001215.4/GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.fna    143726002       25286936
	# Eukaryota/GCA_000001735.2/GCA_000001735.2_TAIR10.1_genomic.fna  119668634       23459830
	# Eukaryota/GCA_000001985.1/GCA_000001985.1_JCVI-PMFA1-2.0_genomic.fna    28643865        3339384
	open my $fh, '<', $file or die "Cannot open $file : $!\n";
	while(my $line=<$fh>) {
		chomp $line;
		my @bits = split /\t/, $line;
		my ($genome_fasta, $length, $n50) = @bits;
		my @genome_fasta_parts = split /\//, $genome_fasta;
		my $assembly = $genome_fasta_parts[1];
		my $fasta_file = $genome_fasta_parts[2];

		#die "assembly = $assembly\n";

		# save
		$data{$assembly}{'length'} = $length;
		$data{$assembly}{'n50'} = $n50;
	}
	return \%data;
}

sub save_coords_info {
	my ($file, $assembly) = @_;
	my %data;
	my %final_features;
	my %gene_to_strand;

	warn "save_coords_info: $file, $assembly\n";

	# AE013599.5	rna-gnl|FlyBase|CG10023-RA	gene-Dmel_CG10023	-	gene_name	ID2	ID3	CDS_range:	19431106-19436596	CDS:	19436548-19436596 19435777-19436479 19435608-19435717 19435314-19435549 19434755-19435250 19434560-19434692 19434370-19434502 19434164-19434308 19433904-19434108 19433703-19433842 19433502-19433641 19432303-19432392 19431792-19432246 19431627-19431737 19431106-19431562
	#AE013599.5	rna-gnl|FlyBase|CG10036-RC	gene-Dmel_CG10036	-	gene_name	ID2	ID3	CDS_range:	20884492-20896492	CDS:	20896475-20896492 20895711-20895866 20886523-20886727 20885202-20885368 20884492-20884761

	open my $fh, '<', $file or die "Cannot open $file : $!\n";
	while(my $line=<$fh>) {
		chomp $line;

		next if($line =~ m/^\n/);
		next if($line eq '');

		my @bits = split /\t/, $line;
		my ($contig, $gene_id, $gene_id2, $strand, $gene_name, $id2, $id3, $range_text, $cds_range, $cds_text, $cds) = @bits;

		# save gene bounderies
		my @cds_range_parts = split /-/, $cds_range;
		my ($start, $stop) = @cds_range_parts;
		$data{$assembly}{$contig}{$start} = $gene_id;
		$data{$assembly}{$contig}{$stop} = $gene_id;

		# save final features
		if(!defined $final_features{$assembly}{$contig}) {
			$final_features{$assembly}{$contig} = $stop
		}
		else {
			if($final_features{$assembly}{$contig} < $stop) {
				$final_features{$assembly}{$contig} = $stop;
			}
		}

		# save strand
		$gene_to_strand{$gene_id} = $strand;
	}
	return (\%data, \%final_features, \%gene_to_strand);	
}

sub save_repeat_info {
	my ($file, $assembly) = @_;
	my %data;
	my %data_names;

	warn "save_repeat_info: $file\n";

	# save repeat_info E.g.
	#     SW   perc perc perc  query       position in query              matching           repeat                position in repeat
	# score   div. del. ins.  sequence    begin    end          (left)   repeat             class/family      begin   end    (left)      ID

	#   769    0.0  0.0  0.0  AE013599.5     16116    16205 (25270731) C rnd-1_family-488   Unknown               (0)    190     101      1  
	#  7404    0.1  0.0  0.1  AE013599.5     16206    17059 (25269877) + rnd-5_family-1189  DNA/TcMar-Pogo          1    853  (1091)      2  
	#   676    0.0  0.0  0.0  AE013599.5     17059    17136 (25269800) + rnd-5_family-1189  DNA/TcMar-Pogo       1511   1588   (356)      3 *
	#  1809    0.0  0.0  0.0  AE013599.5     17130    17355 (25269581) + rnd-5_family-1189  DNA/TcMar-Pogo       1719   1944     (0)      4  
	open my $fh, '<', $file or die "Cannot open $file : $!\n";
	while(my $line=<$fh>) {
		chomp $line;

		next if($line =~ m/^\n/);
		next if($line eq '');
		$line =~ s/^\s+//g;
		
		# ignore header
		next if($line =~ m/^SW\s/);
		next if($line =~ m/^score\s/);

		my @bits = split /\s+/, $line;
		my ($sw_score, $perc_div, $perc_del, $perc_ins, $query_seq, $query_start, $query_stop, $left, $C_or_plus, $repeat, $repeat_class, $repeat_start, $repeat_stop, $id, $non_overlapping) = @bits;

		#An asterisk (*) following the final column (see below example)
  		#indicates that there is a higher-scoring match whose domain partly
		#(<80%) includes the domain of the current match.
		next if($line =~ m/\*$/);

		#die "end here (info = $sw_score, $perc_div, $perc_del, $perc_ins, $query_seq, $query_start, $query_stop, $repeat, $repeat_class, $repeat_start, $repeat_stop, $id, $non_overlapping) : $line\n";

		# save
		$data{$assembly}{$query_seq}{"$query_start-$query_stop"} = $repeat_class;
		$data_names{$assembly}{$query_seq}{"$query_start-$query_stop"} = $repeat;
	}
	return (\%data, \%data_names);
}

sub save_gene_to_sigp4 {
	my ($file, $assembly) = @_;

	warn "save_gene_to_sigp4: $file\n";

	my $line_count = 0;
	my $count = 0;
	my %sec_info;

	open my $fh, '<', $file or die "Cannot open $file: $!\n";
	SIGP: while(my $line=<$fh>) {
		chomp $line;	
		next SIGP if($line =~ m/^#/);
		next SIGP if($line =~ m/^Temporary files in/);

		# split
		my @bits = split /\s+|\t+/, $line;
		my ($name, $Cmax, $pos1, $Ymax, $pos2, $Smax, $pos3, $Smean, $D, $y_or_n, $Dmaxcut, $Networks_used) = @bits;
		my $new_line = join("\t", @bits);
		if($line_count eq 0) { 
			$line_count++;
			warn "The file is split, with 0) $name, 1) $Cmax, 2) $pos1, 4) $Ymax, 5) $pos2, 6) $Smax, 7) $pos3, 8) $Smean, 9) $D, 10) $y_or_n\n"; 
		}

		# Ignore unless secreted
		next SIGP if($y_or_n ne 'Y');

		# Secreted. Count & Print
		$count++;
		
		# Save
		$sec_info{$assembly}{$name} = 1;
	}
	close $fh;
	warn "save_gene_to_sigp4: Signal peptides found: $count\n";
	return \%sec_info;
}
