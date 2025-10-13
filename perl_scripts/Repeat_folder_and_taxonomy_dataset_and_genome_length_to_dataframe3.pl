#!/usr/bin/perl -w
use strict;
use Getopt::Std;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -r <folder of assemblies that have been repeatmodelled, repeatmasked and a summary made> -t <taxonomy_data_for_accessions.tab> -g <genome_lengths.tab>\n
Notes: Makes 2 outfile dataframes, and will probably run an Rscript, once it's written, reading in those outputs\n";
our($opt_r, $opt_t, $opt_g);
getopt('rtg');
die $usage unless ($opt_r && $opt_t && $opt_g);
foreach($opt_r, $opt_t, $opt_g) { die "Cannot open $_ : $!" unless (-e $_); }

# Output files
my $outfile1 = "./RF_12345_Repeat_summary_dataframe1.tab";
my $outfile2 = "./RF_12345_Repeat_summary_dataframe2.tab";
my $outfile3 = "./RF_12345_Repeat_summary_dataframe3.tab";
my $outfile4 = "./RF_12345_Repeat_summary_dataframe4.tab";

# All data
# assembly -> 'kingdom' -> name
# assembly -> species_name -> name
my $assembly_to_taxonomy = &save_taxonomy_data($opt_t);

# Assembly -> length = length
# Assembly -> n50 = n50
my $assembly_metrics = &save_genome_lengths($opt_g);

# Extract just info of interest

my %dataframe_info;
my %kingdom_counts;
my %repeat_type_info; # assembly -> repeat_name -> nt and percent
my %repeat_type_info_summary; # same but condensed (e.g., all DNA transposons joined)

# Save all repeat folders
warn "Save each assembly with repeats in $opt_r...\n";
#my @gff_files = <$opt_r/*>;
my @assemblies_with_repeats_info = <$opt_r/*>;
foreach my $assembly_dir(@assemblies_with_repeats_info) {

	# split e.g. 9.Repeats_all/finished//GCA_000001985.1
	my @assembly_dir_parts = split /\//, $assembly_dir;
	my $assembly = $assembly_dir_parts[-1];
	#warn "assembly = $assembly\n";
	die "No assembly found for $assembly (from folder) in taxonomy file\n" if(!defined $$assembly_to_taxonomy{$assembly});
	die "no kingdom for $assembly (from folder) in taxonomy file\n" if(!defined $$assembly_to_taxonomy{$assembly}{'kingdom'});
	die "no species name for $assembly (from folder) in taxonomy file\n" if(!defined $$assembly_to_taxonomy{$assembly}{'species_name'});
	my $kingdom = $$assembly_to_taxonomy{$assembly}{'kingdom'};
	my $species = $$assembly_to_taxonomy{$assembly}{'species_name'};

	# assembly metrics
	die "no genome length for $assembly\n" if(!defined $$assembly_metrics{$assembly}{'length'});
	die "no n50 for $assembly\n" if(!defined $$assembly_metrics{$assembly}{'n50'});
	my $genome_length = $$assembly_metrics{$assembly}{'length'};
	my $n50 = $$assembly_metrics{$assembly}{'n50'};

	# Find repeat summary made by: for i in */*.out; do perl ~/perl_scripts/RepeatMasker_output_parse.pl -f $i > $i-class_summary.tab; done 
	my @repeat_summary = <$assembly_dir/*-class_summary.tab>;
	die "no repeat summary found in $assembly_dir\n" if(!defined $repeat_summary[0]);
	die "more than 1 repeat summary found in $assembly_dir\n" if(defined $repeat_summary[1]);
	my $repeat_summary_file = $repeat_summary[0];
	my $repeat_summary_info = &save_repeat_summary_info($repeat_summary_file);

	# total repeats
	# $data{'total_repeats'} += $number;
	# individual repeats
	# $data{'repeats'}{$repeat} = $number;
	my $genome_repeat_nt = $$repeat_summary_info{'total_repeats'};
	my $genome_repeat_pc = sprintf("%.2f", (($genome_repeat_nt / $genome_length) *100));
	#die "ive found $repeat_summary_file\n";

	# kingdom count
	$kingdom_counts{$kingdom}++;

	# dataframe info
	$dataframe_info{$assembly}{'kingdom'} = $kingdom;
	$dataframe_info{$assembly}{'species'} = $species;
	$dataframe_info{$assembly}{'genome_length'} = $genome_length;
	$dataframe_info{$assembly}{'n50'} = $n50;
	$dataframe_info{$assembly}{'total_repeat_nt'} = $genome_repeat_nt;
	$dataframe_info{$assembly}{'total_repeat_pc'} = $genome_repeat_pc;

	# dataframe info 2 (Repeat families)
	foreach my $repeat(keys %{$$repeat_summary_info{'repeats'}}) {
		my $tally = $$repeat_summary_info{'repeats'}{$repeat};
		$repeat_type_info{$repeat}{$assembly} = $tally;
	}

	# dataframe info 3 (Repeat families summarised)
	foreach my $repeat(keys %{$$repeat_summary_info{'repeats'}}) {
		my @repeat_parts = split /\//, $repeat;
		my $repeat_type = $repeat_parts[0];
		my $tally = $$repeat_summary_info{'repeats'}{$repeat};
		$repeat_type_info_summary{$repeat_type}{$assembly} += $tally;
	}
}

# warn summary (number of species per kingdom)
foreach my $kingdom(sort keys %kingdom_counts) {
	my $count = $kingdom_counts{$kingdom};
	warn "$kingdom\t$count\n";
}

# open out file handles
open my $ofh1, '>', $outfile1 or die "error: Cannot open $outfile1 : $!";
open my $ofh2, '>', $outfile2 or die "error: Cannot open $outfile2 : $!";
open my $ofh3, '>', $outfile3 or die "error: Cannot open $outfile3 : $!";
open my $ofh4, '>', $outfile4 or die "error: Cannot open $outfile4 : $!";

# make dataframe for plots
print $ofh1 "assembly\tkingdom\tspecies\tgenome_length\tn50\tgenome_repeat_nt\tgenome_repeat_pc\n";
foreach my $assembly(sort keys %dataframe_info) {
	my $kingdom = $dataframe_info{$assembly}{'kingdom'};
	my $species = $dataframe_info{$assembly}{'species'};
	my $genome_length = $dataframe_info{$assembly}{'genome_length'};
	my $n50 = $dataframe_info{$assembly}{'n50'};
	my $total_repeat_nt = $dataframe_info{$assembly}{'total_repeat_nt'};
	my $total_repeat_pc = $dataframe_info{$assembly}{'total_repeat_pc'};

	print $ofh1 "$assembly\t$kingdom\t$species\t$genome_length\t$n50\t$total_repeat_nt\t$total_repeat_pc\n";
}
close $ofh1;

# make dataframe 2 for plots
# header
print $ofh2 "Assembly\tKingdom\tSpecies\tGenome_length";
foreach my $repeat(sort keys %repeat_type_info) {
	print $ofh2 "\t$repeat";
}
print $ofh2 "\n";
# data
foreach my $assembly(sort keys %dataframe_info) {
	my $kingdom = $dataframe_info{$assembly}{'kingdom'};
	my $species = $dataframe_info{$assembly}{'species'};
	my $genome_length = $dataframe_info{$assembly}{'genome_length'};
	print $ofh2 "$assembly\t$kingdom\t$species\t$genome_length";
	foreach my $repeat(sort keys %repeat_type_info) {
		my $tally = 0;
		if(defined $repeat_type_info{$repeat}{$assembly}) { $tally = $repeat_type_info{$repeat}{$assembly}; }
		print $ofh2 "\t$tally";
	}
	print $ofh2 "\n";
}
close $ofh2;

# make dataframe 3 for plots
# header
print $ofh3 "Assembly\tKingdom\tSpecies\tGenome_length";
foreach my $repeat(sort keys %repeat_type_info_summary) {
	print $ofh3 "\t$repeat";
}
print $ofh3 "\n";
# data
foreach my $assembly(sort keys %dataframe_info) {
	my $kingdom = $dataframe_info{$assembly}{'kingdom'};
	my $species = $dataframe_info{$assembly}{'species'};
	my $genome_length = $dataframe_info{$assembly}{'genome_length'};
	print $ofh3 "$assembly\t$kingdom\t$species\t$genome_length";
	foreach my $repeat(sort keys %repeat_type_info_summary) {
		my $tally = 0;
		if(defined $repeat_type_info_summary{$repeat}{$assembly}) { $tally = $repeat_type_info_summary{$repeat}{$assembly}; }
		print $ofh3 "\t$tally";
	}
	print $ofh3 "\n";
}
close $ofh3;

# make dataframe 4 for plots (Repeat type vs kingdom)
print $ofh4 "Repeat_Type\tCount\n";
foreach my $assembly(sort keys %dataframe_info) {
	my $kingdom = $dataframe_info{$assembly}{'kingdom'};
	#my $species = $dataframe_info{$assembly}{'species'};
	#my $genome_length = $dataframe_info{$assembly}{'genome_length'};
	#print $ofh4 "$assembly\t$kingdom\t$species\t$genome_length";
	foreach my $repeat(sort keys %repeat_type_info_summary) {
		my $repeat_type_name = "$repeat\_$kingdom";
		my $tally = 0;
		if(defined $repeat_type_info_summary{$repeat}{$assembly}) { $tally = $repeat_type_info_summary{$repeat}{$assembly}; }
		print $ofh4 "$repeat_type_name\t$tally\n";
	}
}
close $ofh4;

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

sub save_repeat_summary_info {
	my $file = $_[0];
	my %data;

	# save repeat_info E.g.
	# repeat_class	GCA_000001985.1/GCA_000001985.1_JCVI-PMFA1-2.0_genomic.fna.out
	# DNA	21503
	# DNA/Academ	54368
	# DNA/PIF-Harbinger	21000
	# DNA/PiggyBac	33393
	# DNA/TcMar-Fot1	131023
	# DNA/hAT-Tag1	15989
	# LINE/Tad1	39498
	open my $fh, '<', $file or die "Cannot open $file : $!\n";
	while(my $line=<$fh>) {
		chomp $line;
		
		# ignore header
		next if($line =~ m/^repeat_class\t/);

		my @bits = split /\t/, $line;
		my ($repeat, $number) = @bits;

		# total repeats
		$data{'total_repeats'} += $number;

		# individual repeats
		$data{'repeats'}{$repeat} = $number;
	}
	return (\%data);
}
