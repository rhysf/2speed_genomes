package gfffile;
use strict;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);

### rfarrer@broadinstitute.org

sub read_GFF_lines {
	my ($GFF_line, $feature_wanted, $desc_seperator, $desc_column, $desc_replace) = @_;
	my %GFF_info;
	$GFF_info{'next'} = 0;
	if($GFF_line =~ /^#/) {
		$GFF_info{'next'} = 1; 
		$GFF_info{'header'}='Y';
		if($GFF_line =~ m/^\#\#FASTA/) { $GFF_info{'fasta'} = 1; }
		return \%GFF_info;
	}
	if($GFF_line =~ /^$/) { 
		$GFF_info{'next'} = 1; 
		return \%GFF_info;
	}
		
	# initial quality check
	my @bits = split /\t/, $GFF_line;
	if(@bits < 8) {
		warn "Bad GFF with < 8 columns: $GFF_line\n";
		$GFF_info{'next'}=1; 
		return \%GFF_info;
	}

	# Parts continued
	$GFF_info{'contig'}  = $bits[0];
	$GFF_info{'source'}  = $bits[1];
	$GFF_info{'feature'} = $bits[2];
	$GFF_info{'start'}   = $bits[3];
	$GFF_info{'stop'}    = $bits[4];
	$GFF_info{'score'}   = $bits[5];
	$GFF_info{'strand'}  = $bits[6];
	$GFF_info{'frame'}   = $bits[7];
	$GFF_info{'description'} = $bits[8];
	$GFF_info{'CDS'} = "$GFF_info{'start'}-$GFF_info{'stop'}";
	
	# Find feature and parent name
	if(($GFF_info{'feature'} ne $feature_wanted) && ($feature_wanted ne 'all')) { $GFF_info{'feature_wanted'} = 'N'; }
	else {
		$GFF_info{'feature_wanted'} = 'Y';
		my @description_parts = split /$desc_seperator/, $GFF_info{'description'};
		$GFF_info{'feature_parent'} = $description_parts[$desc_column];
		die "read_GFF_lines: Feature parent could not be defined from part $desc_column of $GFF_info{'description'} split by $desc_seperator\n" if(!defined $GFF_info{'feature_parent'});
		$GFF_info{'feature_parent'} =~ s/$desc_replace//g;
	}

	return \%GFF_info;
}

sub gff_to_contig_parent_to_cds_hash {
	my ($input, $feature_wanted, $desc_seperator, $desc_column, $desc_replace) = @_;
	my (%hash_info, %hash_strand);
	warn "gff_to_contig_parent_to_cds_hash: saving $feature_wanted from $input (split col $desc_column by $desc_seperator and remove $desc_replace)...\n";
	open my $fh, '<', $input or die "Cannot open $input: $!\n";
	GFF: while(my $line = <$fh>) {
		chomp $line;
		my $GFF_info = &read_GFF_lines($line, $feature_wanted, $desc_seperator, $desc_column, $desc_replace);
		next GFF if($$GFF_info{'next'} eq 1);
		next GFF if($$GFF_info{'feature_wanted'} eq 'N');
		if(defined $$GFF_info{'fasta'}) {
			warn "FASTA found. Ending subroutine.\n";
			last;
		}
		my ($contig, $feature_parent, $cds, $strand) = ($$GFF_info{'contig'}, $$GFF_info{'feature_parent'}, $$GFF_info{'CDS'}, $$GFF_info{'strand'});

		# Save (avoid leaving a space at the end)
		if(!defined $hash_info{$contig}{$feature_parent}) { $hash_info{$contig}{$feature_parent} = $cds; }
		else { $hash_info{$contig}{$feature_parent} .= " $cds"; }
		$hash_strand{$feature_parent} = $strand;
	}
	close $fh;
	return (\%hash_info, \%hash_strand);
}

sub gff_final_features {
	my $input = $_[0];
	my (%final_feature_per_contig);
	warn "gff_final_features: $input\n";
	open my $fh, '<', $input or die "Cannot open $input: $!\n";
	GFF: while(my $line = <$fh>) {
		chomp $line;
		my @bits = split /\t/, $line;
		next GFF if ((@bits < 8) || ($line =~ m/^\#/));
		my ($contig, $source, $feature, $start, $stop, $score, $strand, $frame, $full_description) = @bits;
		
		# final feature
		if(!defined $final_feature_per_contig{$contig}) { $final_feature_per_contig{$contig} = $stop; }
		elsif($final_feature_per_contig{$contig} < $stop) { $final_feature_per_contig{$contig} = $stop; }
	}
	close $fh;
	return (\%final_feature_per_contig);
}

1;
