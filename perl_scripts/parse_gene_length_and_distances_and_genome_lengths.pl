#!/usr/bin/env perl -w
use strict;
use List::Util qw(sum);
use List::MoreUtils qw(:all);

### r.farrer@exeter.ac.uk

#usage
my $usage = "perl $0 <folder with GCA_ subfolders in>\n";
die $usage unless(@ARGV eq 1);
my $folder = $ARGV[0];

# Load genome sizes
my %genome_sizes;
open my $GSIZES, '<', 'genome_lengths_and_N50.tab' or die "Can't open genome_lengths_and_N50.tab: $!";
while (<$GSIZES>) {
    chomp;
    my ($path, $size, $n50) = split /\t/;
    if ($path =~ /\/(GCA_\d+\.\d+)/) {
        $genome_sizes{$1} = $size;
    }
}
close $GSIZES;

# Open summary file for missing gene length tallies
#open my $MISSING, '>', 'missing_gene_lengths_per_genome.tsv' or die $!;
#print $MISSING "Assembly\tMissingGeneLengthCount\n";

# Open genome-level summary file
open my $SUMMARY, '>', 'genome_level_summary.tsv' or die $!;
print $SUMMARY join("\t", qw(Assembly GenomeSize GeneCount TotalCodingLength CodingFraction MeanGeneLength Median5pFIR Median3pFIR)), "\n";

# Output header
#print join("\t", qw(Assembly Contig GeneID GeneLength UpstreamDistance DownstreamDistance Secreted)), "\n";

# Get all folders
my @folders = glob("$folder/GCA_*");

foreach my $folder (@folders) {
    my ($assembly) = $folder =~ /\/?(GCA_\d+\.\d+)/;
    next unless $assembly;

    my $missing_count = 0;

    # Find the required files using wildcards
    my ($coords_file)    = glob("$folder/*.coords");
    my ($distances_file) = glob("$folder/*distances.tab");
    my ($secreted_file)  = glob("$folder/*.SigP.secreted");

    # Skip if any files are missing
    unless (-e $coords_file && -e $distances_file && -e $secreted_file) {
        warn "Skipping $assembly (missing one or more required files):\n";
        if(! -e $coords_file) { warn "no coords for $folder\n"; }
	if(! -e $distances_file) { warn "no distances for $folder\n"; }
	if(! -e $secreted_file) { warn "no secreted for $folder\n"; }
	next;
    }

    # Read gene lengths from coords
    my %gene_lengths;
    open my $COORDS, '<', $coords_file or die $!;
    while (<$COORDS>) {
        chomp;
        next unless /CDS_range:/;
        my @fields = split /\t/;
        my $gene_id = $fields[1];
	my $cds_block = $fields[10];
	#my ($cds_block) = /CDS:\s+(.+)/ or next;
        my $length = 0;
        foreach my $region (split /\s+/, $cds_block) {
            my ($start, $end) = split /-/, $region;
            $length += abs($end - $start) + 1;
        }
        $gene_lengths{$gene_id} = $length;
    }
    close $COORDS;

    # Read secretion data
    my %secreted;
    open my $SECRETED, '<', $secreted_file or die $!;
    while (<$SECRETED>) {
        chomp;
        my ($gene_id) = split /\t/;
        $secreted{$gene_id} = 1;
    }
    close $SECRETED;

	# Per-genome accumulators
	my (@gene_lengths, @firs_up, @firs_down);
	my $coding_length = 0;
	my $gene_count = 0;

    # Process distances
    open my $DIST, '<', $distances_file or die $!;
    while (<$DIST>) {
        chomp;
        next if /^contig/;  # skip header
        my @fields = split /\t/;
        my ($contig, $gene_id, $up, $down) = @fields[0,1,4,5];

        # Get gene length if it exists
        my $length = exists $gene_lengths{$gene_id} ? $gene_lengths{$gene_id} : 'NA';

        # Count and skip if missing
        if ($length eq 'NA') {
            $missing_count++;
            next;
        }

        # Check if secreted
        my $is_secreted = exists $secreted{$gene_id} ? 1 : 0;

        # Output only valid gene lengths
        $coding_length += $length;
	    push @gene_lengths, $length;
	    push @firs_up, $up;
	    push @firs_down, $down;
	    $gene_count++;
	    #print join("\t", $assembly, $contig, $gene_id, $length, $up, $down, $is_secreted), "\n";
	}
	close $DIST;

	# Write tally
	#print $MISSING "$assembly\t$missing_count\n";

    	my $median_5p = @firs_up ? median(@firs_up) : 0;
	my $median_3p = @firs_down ? median(@firs_down) : 0;
	my $mean_length = $gene_count > 0 ? int($coding_length / $gene_count) : 0;
	my $genome_size = $genome_sizes{$assembly} || 0;
	my $coding_fraction = ($genome_size > 0) ? sprintf("%.4f", $coding_length / $genome_size) : 'NA';

	print $SUMMARY join("\t", $assembly, $genome_size, $gene_count, $coding_length, $coding_fraction, $mean_length, $median_5p, $median_3p), "\n";
}

sub median {
    my @sorted = sort { $a <=> $b } @_;
    my $count = scalar @sorted;
    return $count % 2 ? $sorted[int($count/2)] : ($sorted[$count/2 - 1] + $sorted[$count/2]) / 2;
}
