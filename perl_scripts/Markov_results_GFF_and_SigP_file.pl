#!/usr/bin/perl -w
use strict;
use Getopt::Std;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -m <markov results> -g <gff> -s <SigP4>\n
Optional: -f feature [mRNA]
          -b Feature part after splitting with ; [0]\n";
our($opt_b, $opt_f, $opt_g, $opt_m, $opt_s);
getopt('bfgms');
die $usage unless ($opt_m && $opt_g && $opt_s);
foreach($opt_m, $opt_g, $opt_s) { die "Cannot open $_ : $!" unless (-e $_); }
if(!defined $opt_b) { $opt_b = 0; }
if(!defined $opt_f) { $opt_f = 'mRNA'; }

# Directory
my @gff_file_parts = split /\//, $opt_g;
my $root_dir_gff = $gff_file_parts[0];

# Save markov_coords_for_quadrant_two
# contig -> regions_saved = 1
my $quad2_regions = &save_contig_to_region_saved_for_quad_from_markov($opt_m, 2);
my $quad4_regions = &save_contig_to_region_saved_for_quad_from_markov($opt_m, 4);

# save sigP
my ($secreted, $secreted_first_gene_name) = &save_sigP($opt_s);
my ($secreted_count_quad2, $secreted_count_quad4) = (0, 0);

# Pull out genes from GFF
my $first_gene_print_count = 0;
my $first_GFF_name;
my %genes_in_quad2;
my %genes_in_quad4;
warn "Find regions in GFF...\n";
open my $fh, '<', $opt_g or die "Cannot open $opt_g : $!\n";
while(my $line=<$fh>) {
	chomp $line;
	next if($line =~ m/^$/);
	next if($line =~ m/^#/);
	next if($line =~ m/^\n/);
	my @bits = split /\t/, $line;
	my ($contig, $source, $feature, $start, $stop, $score, $strand, $ignore, $desc) = @bits;

	# lets only save genes
	next if($feature ne $opt_f);

	# get gene
	my @desc_bits = split /\;/, $desc;
	next if(!defined $desc_bits[$opt_b]);
	my $gene = $desc_bits[$opt_b];
	$gene =~ s/^ID=//;
	$gene =~ s/^transcriptId=//;

	if($first_gene_print_count eq 0) {
		#warn "First GFF gene name = $gene\n";
		$first_GFF_name = $gene;
		$first_gene_print_count = 1;
	}

	# check for quad2
	my ($found_quad2) = &check_GFF_line_for_saved_markov_region($quad2_regions, $contig, $start, $stop);
	if($found_quad2 ne 'None') {
		#print "Quad 2\t$contig\t$found_quad2\t$desc\t$gene\n";

		# save genes
		$genes_in_quad2{$gene} = 1;
		if(defined $$secreted{$gene}) { $secreted_count_quad2++; }
	}

	# check for quad4
	my ($found_quad4) = &check_GFF_line_for_saved_markov_region($quad4_regions, $contig, $start, $stop);
	if($found_quad4 ne 'None') {
		#print "Quad 4\t$contig\t$found_quad4\t$desc\t$gene\n";

		# save genes
		$genes_in_quad4{$gene} = 1;
		if(defined $$secreted{$gene}) { $secreted_count_quad4++; }
	}
}
die "ERROR: No genes found\n" if(!defined $first_GFF_name);

# Num genes
my $num_genes_in_quad2 = scalar(keys(%genes_in_quad2));
my $num_genes_in_quad4 = scalar(keys(%genes_in_quad4));
my $pc_quad2_secreted = 0;
my $pc_quad4_secreted = 0;
if($secreted_count_quad2 > 0) { $pc_quad2_secreted = sprintf("%0.2f", (($secreted_count_quad2 / $num_genes_in_quad2) * 100)); }
if($secreted_count_quad4 > 0) { $pc_quad4_secreted = sprintf("%0.2f", (($secreted_count_quad4 / $num_genes_in_quad4) * 100)); }

#print "\n";
#print "GFF\t$opt_g\n";
#print "Markov_results\t$opt_m\n";
#print "SigP\t$opt_s\n";

#print "\n";
#print "Genes_in_Q2\t$num_genes_in_quad2\n";
#print "Secreted_in_Q2\t$secreted_count_quad2\n";
#print "Secreted_in_Q2_pc\t$pc_quad2_secreted\%\n";

#print "\n";
#print "Genes_in_Q4\t$num_genes_in_quad4\n";
#print "Secreted_in_Q4\t$secreted_count_quad4\n";
#print "Secreted_in_Q4_pc\t$pc_quad4_secreted\%\n";
#print "\nFnished\n";

print "Genome\tGFF_name\tSigP_name\tGenes_in_Q2\tSecreted_in_Q2\tSecreted_in_Q2_pc\tGenes_in_Q4\tSecreted_in_Q4\tSecreted_in_Q4_pc\n";
print "$root_dir_gff\t$first_GFF_name\t$secreted_first_gene_name\t$num_genes_in_quad2\t$secreted_count_quad2\t$pc_quad2_secreted\t$num_genes_in_quad4\t$secreted_count_quad4\t$pc_quad4_secreted\n";

sub check_GFF_line_for_saved_markov_region {
	my ($regions_saved, $contig, $start, $stop) = @_;

	# ignore if no contig
	return 'None' if(!defined $$regions_saved{$contig});

	# look for regions
	foreach my $regions(keys %{$$regions_saved{$contig}}) {
		my @regions_parts = split /-/, $regions;
		my ($region_start, $region_stop) = @regions_parts;

		#die "region start = $region_start and region stop = $region_stop\n";

		# check?
		return $regions if(($start >= $region_start) && ($stop <= $region_stop));
	}

	return 'None';
}

sub save_contig_to_region_saved_for_quad_from_markov {
	my ($file, $quadrant_num) = @_;
	my %results;

	warn "save_contig_to_region_saved_for_quad_from_markov: $file ($quadrant_num)...\n";
	open my $fh, '<', $file or die "Cannot open $file : $!\n";
	while(my $line=<$fh>) {
		chomp $line;
		my @bits = split /\t/, $line;
		my ($contig, $region_saved, $quad, $consec_count, $total_count, $prob_quad, $cumulative_prob) = @bits;

		# ignore header
		next if($line =~ m/^Contig\tRegion_saved/);

		# save data for quadrant only
		next if($quad ne $quadrant_num);

		# save each region saved
		my @regions_saved = split / /, $region_saved;
		foreach my $region(@regions_saved) {
			$results{$contig}{$region} = 1;
		}
	}
	close $fh;
	return \%results;
}

sub save_sigP {
	my $file = $_[0];

	# save secreted genes
	my %saved;
	my $first_name;
	my $line_count = 0;
	open my $fh, '<', $file or die "Cannot open $file: $!\n";
	SIGP: while(my $line=<$fh>) {
		chomp $line;	
		next SIGP if($line =~ m/#/);
		my @bits = split /\s+|\t+/, $line;
		#my @bits = split /\t/, $line;
		my ($name, $Cmax, $pos1, $Ymax, $pos2, $Smax, $pos3, $Smean, $D, $y_or_n, $Dmaxcut, $Networks_used) = @bits;
		my $new_line = join("\t", @bits);
		if($line_count eq 0) { 
			$line_count++;
			#print "SigP split, with 0) $name, 1) $Cmax, 2) $pos1, 4) $Ymax, 5) $pos2, 6) $Smax, 7) $pos3, 8) $Smean, 9) $D, 10) $y_or_n\n"; 
			#print "First SigP gene = $name\n";
			$first_name = $name;
		}

		if($y_or_n =~ m/Y/) {
			#if($opt_p eq 'secreted_list') { print "$new_line\n"; }
			$saved{$name} = 1;
		}
	}
	return (\%saved, $first_name);
}
