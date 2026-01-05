#!/usr/bin/env perl -w
use strict;

### r.farrer@exeter.ac.uk

#usage
my $usage = "perl $0 <folder with GCA_ subfolders in> > outfile\n";
die $usage unless(@ARGV eq 1);
my $folder = $ARGV[0];

# Open summary file for missing gene length tallies
open my $MISSING, '>', 'missing_gene_lengths_per_genome.tsv' or die $!;
print $MISSING "Assembly\tMissingGeneLengthCount\n";

# Output header
print join("\t", qw(Assembly Contig GeneID GeneLength UpstreamDistance DownstreamDistance Secreted)), "\n";

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
        print join("\t", $assembly, $contig, $gene_id, $length, $up, $down, $is_secreted), "\n";
    }
    close $DIST;

    # Write tally
    print $MISSING "$assembly\t$missing_count\n";

}
