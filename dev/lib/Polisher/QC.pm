package QC;

use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(fasta_check nt_check);
%EXPORT_TAGS = ( DEFAULT => [qw(&fasta_check)],
                 Both    => [qw(&fasta_check &nt_check)]);


sub fasta_check
{
    my $infile = $_[0];
    if ($infile =~ m/\.gz$/) { ## if a gzip compressed infile
	open(IN,"gunzip -c $infile |") || die "\n\n Cannot open the input file: $infile\n\n";
    }
    else { ## If not gzip compressed
	open(IN,"<$infile") || die "\n\n Cannot open the input file: $infile\n\n";
    }
    my $firstLine = <IN>;
    unless ($firstLine =~ m/^>/) {
	die "\n\n Error: The following file is not in FASTA format:\n $infile\n\n"
    }
    close(IN);
}

sub nt_check
{
    my $infile = $_[0];
    my $bases   = 0;
    my $valid   = 0;
    my $num_seqs = 0;
    if ($infile =~ m/\.gz$/) { ## if a gzip compressed infile
        open(IN,"gunzip -c $infile |") || die "\n\n Cannot open the input file: $infile\n\n";
    }
    else { ## If not gzip compressed
        open(IN,"<$infile") || die "\n\n Cannot open the input file: $infile\n\n";
    }
    while(<IN>) {
	chomp;
	if ($_ =~ m/^>/) {
	    $num_seqs++;
	}
	else {
	    my $sequence = $_;
	    $bases += length($sequence);
	    $valid += $sequence =~ tr/ATGCatgc/ATGCatgc/;
	}
    }
    close(IN);
    my $percent = ($valid / $bases) * 100;
    unless ( $percent > 98 ) {
	die "\n\n Error: This file does not appear to contain nucleotides as only $percent % of characters are A,T,G, or C\n\n";
    }
    return($num_seqs);
}

sub thread_check
{
    my $seqs    = $_[0];
    my $threads = $_[1];
    if ($seqs >= $threads) {
	return ($threads);
    }
    else {
	print " @@ Warning @@\n Thread count reduced from $threads to $seqs becasue\n you do not have enough sequences to parallelize by $threads threads\n\n";
	return ($seqs);
    }
}

1;
