#!/usr/bin/perl -w

use strict;
use FindBin;
use Cwd 'abs_path';
use lib abs_path("$FindBin::Bin/../lib");
use File::Basename;
use Polisher::QC qw(:Both);
use Polisher::Alignment qw(:Both);
use Polisher::Format;

my $reference_offset = $ARGV[0];
my $reference_file   = $ARGV[1];
my $read_file        = $ARGV[2];

my ($reference,$read);
my (%three_frame,%positive_scoring_regions,%framer,%tie);

## Needleman-Wunsch Scoring Parameters
my $MATCH    =  1;
my $MISMATCH = -1;
my $GAP      = -1;

open(IN,"<$reference_file") || die "\n Cannot open the reference file: $reference_file\n";
while(<IN>) {
    chomp;
    $reference = $reference . $_;
}
close(IN);

open(IN,"<$read_file") || die "\n Cannot open the read file: $read_file\n";
while(<IN>) {
    chomp;
    $read = $read . $_;
}
close(IN);

## Trim down the reference using the reference_offset
my $trmd_reference = substr $reference, ($reference_offset-1);
my ($ref_min,$ref_max,$q_begin_pos,$q_end_pos) = (10**9,0,10**9,0);

## NW alignment for each frame
for (my $frame=1; $frame <=3; $frame++) {
    my ($q_pos,$r_pos) = (1,1);
    my $query_peptide = Alignment::translate($read, $frame);
    $three_frame{$frame} = $query_peptide;
    my ($q_align,$aln_string,$r_align,$score) = Alignment::needleman_wunsch($query_peptide, $trmd_reference, $MATCH, $MISMATCH, $GAP);
    for (my $i=0; $i <= length($q_align)-3; $i++) {
	my $q_window = substr $q_align, $i, 3;
	my $r_window = substr $r_align, $i, 3;
	my $score = Alignment::score($q_window, $r_window, $MATCH, $MISMATCH, $GAP);
	if ($q_pos == 1 && $r_pos == 1) {
	    my $q_anchor = substr $q_window, 0, 1;
	    my $r_anchor = substr $r_window, 0, 1;
	    if ($score > 0) {
		$positive_scoring_regions{$frame}{$r_pos} = $q_pos;
		if ($r_pos < $ref_min) { $ref_min = $r_pos; $q_begin_pos = $q_pos; }
		if ($r_pos > $ref_max) { $ref_max = $r_pos; $q_end_pos   = $q_pos; }
	    }
	    if ($q_anchor ne "-") { $q_pos++; }
            if ($r_anchor ne "-") { $r_pos++; }
	    $q_anchor = substr $q_window, 1, 1;
            $r_anchor = substr $r_window, 1, 1;
            if ($score > 0) {
                $positive_scoring_regions{$frame}{$r_pos} = $q_pos;
		if ($r_pos < $ref_min) { $ref_min = $r_pos; $q_begin_pos = $q_pos; }
		if ($r_pos > $ref_max) { $ref_max = $r_pos; $q_end_pos   = $q_pos; }
	    }
            if ($q_anchor ne "-") { $q_pos++; }
            if ($r_anchor ne "-") { $r_pos++; }
	}
	else {
	    my $q_anchor = substr $q_window, 1, 1;
            my $r_anchor = substr $r_window, 1, 1;
            if ($score > 0) {
                $positive_scoring_regions{$frame}{$r_pos} = $q_pos;
		if ($r_pos < $ref_min) { $ref_min = $r_pos; $q_begin_pos = $q_pos; }
		if ($r_pos > $ref_max) { $ref_max = $r_pos; $q_end_pos   = $q_pos; }
	    }
            if ($q_anchor ne "-") { $q_pos++; }
            if ($r_anchor ne "-") { $r_pos++; }
	}
    }
}

foreach my $i (sort {$a<=>$b} keys %positive_scoring_regions) {
    foreach my $j (sort {$a<=>$b} keys %{$positive_scoring_regions{$i}}) {
	if (exists $framer{$positive_scoring_regions{$i}{$j}}) {
	    if (exists $tie{$framer{$positive_scoring_regions{$i}{$j}}}) {
		$tie{$positive_scoring_regions{$i}{$j}} = $tie{$framer{$positive_scoring_regions{$i}{$j}}} . "_" . $i;
	    }
	    else {
		$tie{$positive_scoring_regions{$i}{$j}} = $framer{$positive_scoring_regions{$i}{$j}} . "_" . $i;
	    }
	}
	else {
	    $framer{$positive_scoring_regions{$i}{$j}} = $i;
	}
    }
}

for (my $pos = $q_begin_pos; $pos <= $q_end_pos; $pos++) {
    if (exists $framer{$pos}) {
	my $aa = substr $three_frame{$framer{$pos}}, $pos-1, 1;
	if (exists $tie{$pos}) {
	    my ($previous_frame,$next_frame);
	    for (my $i=1;$i<$pos-$q_begin_pos;$i++) {
		if (exists $framer{$pos-$i}) {
		    $previous_frame = $framer{$pos-$i};
		    last;
		}
	    }
	    for(my $i=1;$i<$q_end_pos - $pos;$i++) {
		if (exists $framer{$pos+$i}) {
                    $next_frame = $framer{$pos-$i};
                    last;
                }
            }
	    if (defined $previous_frame) { print "$previous_frame"; }
	    print " <-> ";
	    if (defined $next_frame) { print "$next_frame"; }
	    print "\n";
	}
	else {
	    # print "$aa ";
	}
    }
    else {
	# print "? ";
    }
}

print "\n";

exit 0;
