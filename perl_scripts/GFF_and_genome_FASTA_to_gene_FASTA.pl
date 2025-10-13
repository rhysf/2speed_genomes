#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/perl_modules";
use MutationTools::read_FASTA;
use MutationTools::read_GFF;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "perl $0 -f <genome.fasta> -g <GFF> > out.fasta
Optional: -a\tFeature in GFF for sequence (exon/CDS/gene) [CDS]
          -s\tSeperator in description for gene names (\" ; etc) [;]
          -c\tColumn that says which parent/gene it belongs to the description [0]
          -r\tRemove additional comments in column [Parent=]
          -b\tFeature in GFF for all the ids (exon/CDS/gene) [mRNA]\n";
our($opt_a, $opt_b, $opt_c, $opt_f, $opt_g, $opt_s, $opt_r, $opt_u);
getopt('abcfgsru');
die $usage unless ($opt_f && $opt_g);
if(!defined $opt_a) { $opt_a = 'CDS'; }
if(!defined $opt_b) { $opt_b = 'mRNA'; }
if(!defined $opt_s) { $opt_s = ';'; }
if(!defined $opt_c) { $opt_c = 0; }
if(!defined $opt_u) { $opt_u = 0; }
if(!defined $opt_r) { $opt_r = 'Parent='; }

# Save sequences
my ($sequences, $descriptions, $order) = fastafile::fasta_id_to_seq_hash($opt_f);

# Save GFF
my ($gff, $strand) = gfffile::gff_to_parent_to_cds_hash($opt_g, $opt_a, $opt_s, $opt_c, $opt_r);
my ($gff_to_contig, $strand_ignore) = gfffile::gff_to_contig_parent_to_cds_hash($opt_g, $opt_a, $opt_s, $opt_c, $opt_r);

# Save GFF2 should be formatted like so: ID=7000010424395909;Alias=CNBG_4576T0;Parent=7000010424395908;Name=peroxisomal copper amine oxidase
my ($gff2, $strand2) = gfffile::gff_to_parent_to_cds_hash($opt_g, $opt_b, '@@@@', 0, '@@@@');

# Find GFF in GFF2
my %gff_to_description;
foreach my $full_description(keys %{$gff2}) {
	my @desc_parts = split /$opt_s/, $full_description;
	my ($id, $alias, $parent, $name) = @desc_parts;
	die "No $id found in GFF\n" if(!defined $id);
	if(!defined $alias) { $alias = ''; }
	if(!defined $parent) { $parent = ''; }
	if(!defined $name) { $name = ''; }
	$id =~ s/^ID=//;
	$alias =~ s/^Alias=//;
	$alias =~ s/T[0123456789]+$//;
	$parent =~ s/^Parent=//;
	$name =~ s/^Name=/name=\"/;
	$name .= '"';
	my $description_for_fasta = "$id $parent $alias $name";
	if(defined $$gff{$id}) { $gff_to_description{$id} = $description_for_fasta; }
	else {
		warn "Did not find $id from $full_description in the GFF ($opt_b) also in $opt_a\n"; 
		last;
	}
}

# Pull out all the exons
my %CDS_extracted;
warn "Extracting all CDS from $opt_f...\n";
CONTIG: foreach my $contig(sort keys %{$sequences}) {
	GENES: foreach my $gene_id(sort keys %{$$gff_to_contig{$contig}}) {
		my $gene_description = $gene_id;
		if(defined $gff_to_description{$gene_id}) { $gene_description = $gff_to_description{$gene_id}; }
		#die "gene description not defined $contig -> $gene_id\n" if(!defined $gene_description);
		my $start_stop = $$gff_to_contig{$contig}{$gene_id};

		# new method
		#warn "gene = $gene_id\n";
		$CDS_extracted{$gene_description} = &extract_gene_from_coords_new_negative($$sequences{$contig}, $start_stop, $$strand{$gene_id}, $opt_u, $contig); 

		# delete failed attempts
		if($CDS_extracted{$gene_description} eq 'N') { delete $CDS_extracted{$gene_description}; }
	}
}

# Print
my $print = fastafile::fasta_hash_print_simple(\%CDS_extracted);

sub order_exons_and_star_stop_within_exons {
	my ($exon_string, $strand) = @_;
	my @exons = split /\s/, $exon_string;

	# Make sure each exon is from small number to big number
	for(my $h=0; $h < scalar(@exons); $h++) {
		my @terminals = split /-/, $exons[$h];
		if($terminals[1] < $terminals[0]) {
			warn "Unusual order in exons found: @exons ... check this looks fine then delete comment: $exons[$h]\n";
			my $terminal_temp = $terminals[0];
			$terminals[0] = $terminals[1];
			$terminals[1] = $terminal_temp;
			$exons[$h] = ($terminals[0] . '-' . $terminals[1]);
			warn "new order: $exons[$h]\n";
		}
	}

	# Order exons from end (big) to start (small)
	if($strand eq '-') {
		#warn "strand = -\n";
		#warn "exons: @exons\n";
		#next;
		# Alphanumeric for [\d+-\d+\s]+
		#my @exons_negative = sort { $b <=> $a || $b cmp $a } @exons;
		#my @exons_negative = sort { $b cmp $a } @exons;
		my @exons_negative = sort { $b cmp $a } @exons;
		@exons = @exons_negative;
		#warn "to @exons\n";
	}
	return @exons;
}

sub extract_gene_from_coords_new_negative {
	my ($seq, $exon_string, $strand, $number_nt_to_add_to_terminal_exons, $contig) = @_;
	my ($extracted_gene);

	my $length_of_contig = length($seq);

	# Order exons (- strand goes from small to big for substr purposes) and order small to big start/stop within exons 
	#my @exons = split /\s/, $exon_string;
	my @exons = &order_exons_and_star_stop_within_exons($exon_string, $strand);

	#next if($strand eq '-');
	#warn "gene with $strand strand\n";
	#foreach my $exon(@exons) {
	EXONS: for(my $h=0; $h < scalar(@exons); $h++) {
		my @terminals = split /-/, $exons[$h];
		#warn "exons/terminals: $h) @terminals from all of them: @exons\n";

		# Terminal exons add (e.g. fake UTRs)
		if($number_nt_to_add_to_terminal_exons ne 0) {
			if($h eq 0) {
				#warn "Changing first one: $terminals[0]-$terminals[1]\n";
				if($strand eq '+') {
					$terminals[0] -= $number_nt_to_add_to_terminal_exons;
					if($terminals[0] < 1) { $terminals[0] = 1; }
				} 
				# - strand needs stop position changed e.g. 5369-5590
				else {
					$terminals[1] += $number_nt_to_add_to_terminal_exons;
					if($terminals[1] > $length_of_contig) { $terminals[1] = $length_of_contig; }
				}
				#warn "to $terminals[0]-$terminals[1]\n";
			}

			if($h eq (scalar(@exons) - 1)) {
				my $last_exon_start = (scalar(@terminals) - 2); 
				my $last_exon_stop = (scalar(@terminals) - 1); # start/stop might be inaccurate names depending on strand
				#warn "changing last one: $terminals[$last_exon_start]-$terminals[$last_exon_stop]\n";
				if($strand eq '+') {
					$terminals[$last_exon_stop] += $number_nt_to_add_to_terminal_exons;
					if($terminals[$last_exon_stop] > $length_of_contig) { $terminals[$last_exon_stop] = $length_of_contig; }
				} 
				# - strand needs stop position changed e.g. 5369-5590
				else {
					$terminals[$last_exon_start] -= $number_nt_to_add_to_terminal_exons;
					if($terminals[$last_exon_start] < 1) { $terminals[$last_exon_start] = 1; }
				}
				#warn "to $terminals[$last_exon_start]-$terminals[$last_exon_stop]\n";
			}
		}
		#warn "after adding fake UTRs, my terminals are @terminals\n";

		# Extract sequence
		for(my $i=0; $i< scalar(@terminals); $i+=2) { 

			# +1 to turn array position to genome position
			my $pull_start = ($terminals[$i] - 1);
			my $pull_end = ($terminals[($i + 1)] - $terminals[$i] + 1);
			#warn "for @terminals: I want to pull out seq from $pull_start to $pull_end\n";

			# check sequence goes up to this?
			if($length_of_contig < ($pull_start + $pull_end)) {
				warn "WARNING: exons $pull_start-$pull_end (strand $strand) are outside of contig $contig length ($length_of_contig). Therefore ignoring this gene.\n";
				$extracted_gene = 'N';
				last EXONS;
			}

			my $exon_sequence = substr $seq, $pull_start, $pull_end;
			if($strand eq '-') { $exon_sequence = fastafile::reverse_compliment($exon_sequence); }
			$extracted_gene .= $exon_sequence;
		}
	}
	#die "end here\n";
	# New. Uppercase
	$extracted_gene = uc $extracted_gene;
	#die "end here with $extracted_gene\n";
	return $extracted_gene;
}
