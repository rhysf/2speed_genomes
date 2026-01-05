#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/perl_modules";
use MutationTools::read_GFF;
use Data::Dumper;
use LWP::Simple;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -n <assembly_summary_genbank.txt> -g <folder of GFFs downloaded from NCBI>\n";
our($opt_n, $opt_g);
getopt('ng');
die $usage unless ($opt_n && $opt_g);
die "Cannot open $opt_n : $!" unless (-e $opt_n);

# Output directories
`mkdir Archaea`;
`mkdir Bacteria`;
`mkdir Eukaryota`;
`mkdir Viruses`;
`mkdir Failed_QC`;
`mkdir Unknown_superkingdom`;

# Save accession to TaxIDs
my $accession_to_taxid = &save_taxID_from_assembly_summary_genbank($opt_n);

# Go through all GFF files
warn "Save each Gff file in $opt_g (suffix .gff)...\n";
my @gff_files = <$opt_g/*.gff>;

# Loop through all gff files
# 1. check there are Genes in GFF. If none warn and ignore, finding taxid, and saving a copy to the euk, or bacterial folder
GFF: foreach my $gff(@gff_files) {

	# How many genes?
	warn "opening $gff...\n";
	#my $feature_wanted = 'gene';
	#my $desc_seperator = ';';
	#my $desc_column = 0;
	#my $desc_replace = 'ID=';
	#my ($gene_info, $gene_strand) = gfffile::gff_to_parent_to_cds_hash($gff, $feature_wanted, $desc_seperator, $desc_column, $desc_replace); 
	#my $gene_count = scalar(keys(%{$gene_info}));
	#warn "\t$gene_count genes found\n";

	# QC 1 = move to Failed QC if too low!
	# less than 182 too low! (https://www.nature.com/news/2006/061009/full/news061009-10.html)
	#if($gene_count < 182) {
	#	warn "Failed QC: $gff\n";
	#	`mv $gff ./Failed_QC/`;
	#	next GFF;
	#}

	# Next look up TaxID (e.g. .//GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.gff)
	my @gff_parts = split /\_/, $gff;
	$gff_parts[0] =~ s/\.\/+//;
	my $accession = ($gff_parts[0] . '_' . $gff_parts[1]);
	die "No accession from $opt_n found for gff $gff with accession $accession\n" if(!defined $$accession_to_taxid{$accession});
	my $taxID = $$accession_to_taxid{$accession};
	warn "tax id = $taxID\n";

	# Find superkingdom
	my $url = 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=' . $taxID;
	my $superkingdom = &extract_superkingdom_from_NCBI_taxonomy_HTML($url);
	#warn "superkingdom = $superkingdom\n";

	# move it to suitable folder for future work (gene count summaries, maybe SigP?, something to do with 5' and 3' distance) ! (can't compare bacteria to eukaryotes)
	die "What superkingdom has been found for $gff with taxid $taxID: $superkingdom\n" if($superkingdom !~ m/Archaea|Bacteria|Eukaryota|Viruses|Unknown_superkingdom/);
	`mv $gff ./$superkingdom/`;
	#die "end here\n";
}

sub extract_superkingdom_from_NCBI_taxonomy_HTML {
	my $url = $_[0];
	my $superkingdom = 'Unknown_superkingdom';

	# Get HTML
	my $HTML = get($url);

	# Parse for superkingdom
	my @start = split /TITLE="superkingdom"/, $HTML;
	$start[1] =~ s/^>//; # if thats the end of the HTML AHREF
	$start[1] =~ s/^ ALT="superkingdom">//; # or it might have that
	my @end = split /<\/[aA]>;/, $start[1];
	
	if(defined $end[0]) { $superkingdom = $end[0]; }
	#print "$HTML\n";
	return $superkingdom;
}

sub save_taxID_from_assembly_summary_genbank {
	my $file = $_[0];
	my %assembly_assession_to_TAXID;
	#warn "save_taxID_from_assembly_summary_genbank: $file\n";
	open my $fh, '<', $file or die "Cannot open $file : $!\n";
	while(my $line=<$fh>) {
		chomp $line;
		my @bits = split /\t/, $line;
		#   See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt for a description of the columns in this file.
		#Column  6: "taxid"
		#Taxonomy ID: the NCBI taxonomy identifier for the organism from which the
		#genome assembly was derived. The NCBI Taxonomy Database is a curated
		#classification and nomenclature for all of the organisms in the public
		#sequence databases. The taxonomy record can be retrieved from the NCBI
		#Taxonomy resource:
		#https://www.ncbi.nlm.nih.gov/taxonomy/

		#Column  7: "species_taxid"
		#Species taxonomy ID: the NCBI taxonomy identifier for the species from which
		#the genome assembly was derived. The species taxid will differ from the
		#organism taxid (column 6) only when the organism was reported at a sub-
		#species or strain level.
		my ($assembly_accession, $bioproject, $biosample, $wgs_master, $refseq_category, $taxid, $species_taxid, $organism_name) = @bits;
		$assembly_assession_to_TAXID{$assembly_accession} = $taxid;
	}
	return \%assembly_assession_to_TAXID;
}
