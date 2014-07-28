package Alignment;

use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(needleman_wunsch);
%EXPORT_TAGS = ( DEFAULT => [qw(&needleman_wunsch)],
		 Both    => [qw(&needleman_wunsch &score)]);

sub needleman_wunsch
{
    my $seq1 = $_[0];
    my $seq2 = $_[1];
    my $MATCH    = $_[2];
    my $MISMATCH = $_[3];
    my $GAP      = $_[4];
    
    # initialization
    my @matrix;
    $matrix[0][0]{score}   = 0;
    $matrix[0][0]{pointer} = "none";
    for(my $j = 1; $j <= length($seq1); $j++) {
	$matrix[0][$j]{score}   = $GAP * $j;
	$matrix[0][$j]{pointer} = "left";
    }
    for (my $i = 1; $i <= length($seq2); $i++) {
	$matrix[$i][0]{score}   = $GAP * $i;
	$matrix[$i][0]{pointer} = "up";
    }

    # fill
    for(my $i = 1; $i <= length($seq2); $i++) {
	for(my $j = 1; $j <= length($seq1); $j++) {
	    my ($diagonal_score, $left_score, $up_score);

        # calculate match score
	    my $letter1 = substr($seq1, $j-1, 1);
	    my $letter2 = substr($seq2, $i-1, 1);                            
	    if ($letter1 eq $letter2) {
		$diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
	    }
	    else {
		$diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
	    }

        # calculate gap scores
	    $up_score   = $matrix[$i-1][$j]{score} + $GAP;
	    $left_score = $matrix[$i][$j-1]{score} + $GAP;

        # choose best score
	    if ($diagonal_score >= $up_score) {
		if ($diagonal_score >= $left_score) {
		    $matrix[$i][$j]{score}   = $diagonal_score;
		    $matrix[$i][$j]{pointer} = "diagonal";
		}
		else {
		    $matrix[$i][$j]{score}   = $left_score;
		    $matrix[$i][$j]{pointer} = "left";
		}
	    } else {
		if ($up_score >= $left_score) {
		    $matrix[$i][$j]{score}   = $up_score;
		    $matrix[$i][$j]{pointer} = "up";
		}
		else {
		    $matrix[$i][$j]{score}   = $left_score;
		    $matrix[$i][$j]{pointer} = "left";
		}
	    }
	}
    }

    # trace-back

    my $align1 = "";
    my $align2 = "";

    # start at last cell of matrix
    my $j = length($seq1);
    my $i = length($seq2);
    
    while (1) {
	last if $matrix[$i][$j]{pointer} eq "none"; # ends at first cell of matrix

	if ($matrix[$i][$j]{pointer} eq "diagonal") {
	    $align1 .= substr($seq1, $j-1, 1);
	    $align2 .= substr($seq2, $i-1, 1);
	    $i--;
	    $j--;
	}
	elsif ($matrix[$i][$j]{pointer} eq "left") {
	    $align1 .= substr($seq1, $j-1, 1);
	    $align2 .= "-";
	    $j--;
	}
	elsif ($matrix[$i][$j]{pointer} eq "up") {
	    $align1 .= "-";
	    $align2 .= substr($seq2, $i-1, 1);
	    $i--;
	}    
    }

    $align1 = reverse $align1;
    $align2 = reverse $align2;
    my $aln_string = "";
    my @A1 = split(//, $align1);
    my @A2 = split(//, $align2);
    for (my $i=0; $i<length($align1); $i++) {
	if ($A1[$i] eq $A2[$i]) {
	    $aln_string = $aln_string . $A1[$i];
	}
	else {
	    $aln_string = $aln_string . " ";
	}
    }
    my $score = score($align1,$align2,$MATCH,$MISMATCH,$GAP);
    return($align1,$aln_string,$align2,$score);
}

sub score
{
    my $align1   = $_[0];
    my $align2   = $_[1];
    my $MATCH    = $_[2];
    my $MISMATCH = $_[3];
    my $GAP      = $_[4];
    my $score = 0;
    my @A1 = split(//, $align1);
    my @A2 = split(//, $align2);
    for (my $i=0; $i<length($align1); $i++) {
        if ($A1[$i] eq $A2[$i]) {
	    $score += $MATCH;
        }
        else {
            if ($A1[$i] eq "-" || $A2[$i] eq "-") {
		$score += $GAP;
	    }
	    else {
		$score += $MISMATCH;
	    }
	}
    }
    return($score);
}

sub translate
{
    my %codon = (
    'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'TTA' => 'L',
    'TTG' => 'L', 'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
    'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R', 'AGA' => 'R',
    'AGG' => 'R', 'AAA' => 'K', 'AAG' => 'K', 'AAT' => 'N', 'AAC' => 'N',
    'ATG' => 'M', 'GAT' => 'D', 'GAC' => 'D', 'TTT' => 'F', 'TTC' => 'F',
    'TGT' => 'C', 'TGC' => 'C', 'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P',
    'CCG' => 'P', 'CAA' => 'Q', 'CAG' => 'Q', 'TCT' => 'S', 'TCC' => 'S',
    'TCA' => 'S', 'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S', 'GAA' => 'E',
    'GAG' => 'E', 'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
    'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G', 'TGG' => 'W',
    'CAT' => 'H', 'CAC' => 'H', 'TAT' => 'Y', 'TAC' => 'Y', 'ATT' => 'I',
    'ATC' => 'I', 'ATA' => 'I', 'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V',
    'GTG' => 'V', 'TAG' => '*', 'TAA' => '*', 'TGA' => '*' );
    my $seq   = $_[0];
    my $frame = $_[1];
    $seq = uc $seq;
     my $peptide = "";
    if ($frame > 3) { $seq = revcomp($seq);}
    if ($frame == 2 || $frame == 5) {
        $seq = substr $seq, 1;
    }
    elsif ($frame == 3 || $frame == 6) {
        $seq = substr $seq, 2;
    }
    my $size = length($seq);
   for (my $i=0; $i<$size-3; $i+=3) {
        my $trip = substr $seq, 0, 3;
	$peptide = $peptide . $codon{$trip};
	if (length($seq) >= 3) {
	    $seq = substr $seq, 3;
        }
    }
    return($peptide);
}

sub revcomp
{
    my $seq = $_[0];
    my $rev_comp = scalar reverse $seq;
    $rev_comp =~ tr/ATGCatgc/TACGTACG/;
    return($rev_comp);
}

1;
