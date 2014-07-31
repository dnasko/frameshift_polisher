#!/usr/bin/perl -w

use strict;
use FindBin;
use Getopt::Long;
use File::Basename;
use Pod::Usage;
use Cwd 'abs_path';
use lib abs_path("$FindBin::Bin/../lib");
use File::Basename;
use Polisher::QC qw(:Both);
use Polisher::Alignment qw(:Both);
use Polisher::Format;

my ($reference_offset,$reference,$read,$sense,$name);

GetOptions (
    "o|offset=s"=>\$reference_offset,
    "r|ref=s"=>\$reference,
    "q|query=s"=>\$read,
    "s|sense=s"=>\$sense,
    "n|name=s"=>\$name);

# my ($reference,$read);
my ($frame_start,$frame_end);
my (%three_frame,%positive_scoring_regions,%pos_scr_rgns_inv,%framer,%tie);

## Needleman-Wunsch Scoring Parameters
my $MATCH    =  1;
my $MISMATCH = -1;
my $GAP      = -1;

# open(IN,"<$reference_file") || die "\n Cannot open the reference file: $reference_file\n";
# while(<IN>) {
#     chomp;
#     $reference = $reference . $_;
# }
# close(IN);

# open(IN,"<$read_file") || die "\n Cannot open the read file: $read_file\n";
# while(<IN>) {
#     chomp;
#     $read = $read . $_;
# }
# close(IN);

if ($sense == 0) {
    $frame_start = 1;
    $frame_end   = 3;
}
elsif ($sense == 1) {
    $frame_start = 4;
    $frame_end   = 6;
}
else {
    die "\n Error: frame information missing!\n";
}

## Trim down the reference using the reference_offset
my $trmd_reference = substr $reference, ($reference_offset-1);
my ($ref_min,$ref_max,$q_begin_pos,$q_end_pos) = (10**9,0,10**9,0);

## NW alignment for each frame
for (my $frame=$frame_start; $frame <=$frame_end; $frame++) {
    my ($q_pos,$r_pos) = (1,1);
    my $query_peptide = Alignment::translate($read, $frame);
    $three_frame{$frame} = $query_peptide;
    my ($q_align,$aln_string,$r_align,$score) = Alignment::needleman_wunsch($query_peptide, $trmd_reference, $MATCH, $MISMATCH, $GAP);
    ## Slide through the NW alignment
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

## Getting our arms around the ties . . .
foreach my $i (sort {$a<=>$b} keys %positive_scoring_regions) {  ## foreach frame
    foreach my $j (sort {$a<=>$b} keys %{$positive_scoring_regions{$i}}) {
	if (exists $framer{$positive_scoring_regions{$i}{$j}}) {
	    if (exists $tie{$positive_scoring_regions{$i}{$j}}) {
		$tie{$positive_scoring_regions{$i}{$j}} = $tie{$framer{$positive_scoring_regions{$i}{$j}}} . "_" . $i;
	    }
	    else {
		$tie{$positive_scoring_regions{$i}{$j}} = $framer{$positive_scoring_regions{$i}{$j}} . "_" . $i;
	    }
	}
	else {
	    $framer{$positive_scoring_regions{$i}{$j}} = $i;
	}
	$pos_scr_rgns_inv{$i}{$positive_scoring_regions{$i}{$j}} = $j; 
    }
}

print STDOUT "TIES:\n";
foreach my $i (sort {$a<=>$b} keys %tie) {
    # my @frames = split(/_/, $tie{$i});
    # foreach my $f (@frames) {
    # 	if (exists $pos_scr_rgns_inv{$f}{$i}) {
    # 	    print STDOUT "$i\t$pos_scr_rgns_inv{$f}{$i}\t$f\n";
    # 	}
    # 	else {
    # 	    print STDOUT "MISSING: $i\t$f\n";
    # 	}
    # }
    print "$i\t$tie{$i}\n";
}

## Running thorugh %framer to make sure each position has a call
for (my $pos = $q_begin_pos; $pos <= $q_end_pos; $pos++) {
    if (exists $framer{$pos}) {
	my $aa = substr $three_frame{$framer{$pos}}, $pos-1, 1;
	if (exists $tie{$pos}) {
	    my ($previous_frame,$next_frame) = find_prev_next($pos,$q_begin_pos);
	    if (defined $previous_frame && defined $next_frame) {
		if ($previous_frame == $next_frame) {
		    $framer{$pos} = $previous_frame;
		}
		else {
		    my $prev_ref_pos = $pos_scr_rgns_inv{$framer{$pos-1}}{$pos-1};
		    my @ties = split(/_/, $tie{$pos});
		    my $found_flag = 0;
		    foreach my $tie_frame (@ties) {
			print STDERR "TIE: $tie_frame\nPOS: $pos\n";
			my $proposed_ref_pos = $pos_scr_rgns_inv{$tie_frame}{$pos};
			if ($prev_ref_pos+1 == $proposed_ref_pos) {
			    $framer{$pos} = $tie_frame;
			    $found_flag = 1;
			    last;
			}
		    }
		    if ($found_flag == 0) {
			die "\n @@@@@@@ Well that didn't work!@@@@@@@\n\n";
		    }
		}
	    }
	    elsif (defined $previous_frame) {
		my @ties = split(/_/, $tie{$pos});
		my %t;
		foreach my $tie_frame (@ties) {
		    $t{$tie_frame} = 1;
		}
		if (exists $t{$previous_frame}) {
		    $framer{$pos} = $previous_frame;
		}
		else {
		    $framer{$pos} = $ties[0];
		}
	    }
	    elsif (defined $next_frame) {
		my @ties = split(/_/, $tie{$pos});
		my %t;
		foreach my $tie_frame (@ties) {
                    $t{$tie_frame} = 1;
		}
		if (exists $t{$next_frame}) {
                    $framer{$pos} = $next_frame;
		}
		else {
                    $framer{$pos} = $ties[0];
                }
            }
	    else {
		my @ties = split(/_/, $tie{$pos});
		$framer{$pos} = $next_frame;
	    }
	}
    }
    else {
	my ($previous_frame,$next_frame) = find_prev_next($pos,$q_begin_pos);
	if (defined $previous_frame && defined $next_frame) {
	    if ($previous_frame == $next_frame) {
		$framer{$pos} = $previous_frame;
	    }
	    else {
		$framer{$pos} = $next_frame;
	    }
	}
	elsif (defined $previous_frame) {
	    $framer{$pos} = $previous_frame;
	}
	elsif (defined $next_frame) {
	    $framer{$pos} = $next_frame;
	}
	else {
	    ## Hold on, it's probably just the first position...
	    # die "Can't find a previous or upcoming frame:\nPOS:$pos\nBEGIN: $q_begin_pos\nEND:   $q_end_pos\n$read\n\n";
	}
    }
}

my ($corrected_peptide,$frame_string);
for (my $pos = $q_begin_pos; $pos <= $q_end_pos; $pos++) {
    if (exists $framer{$pos}) {
	my $aa;
	if (defined $framer{$pos}) {
	    $aa = substr $three_frame{$framer{$pos}}, $pos-1, 1;
	}
	else {
	    $aa = substr $three_frame{$frame_start}, $pos-1, 1;
	}
# if ($aa eq "*") {
	#     if ($framer{$pos-1} == $framer{$pos} && $framer{$pos} == $framer{$pos+1}) {
	# 	$framer{$pos} = $framer{$pos-1};
		
	#     }
	#     elsif ($framer{$pos} != $framer{$pos+1}) {
	# 	$framer{$pos} = $framer{$pos+1};
	#     }
	#     elsif ($framer{$pos} != $framer{$pos-1}) {
	# 	$framer{$pos} =$framer{$pos-1};
	#     }
	#     $aa = substr $three_frame{$framer{$pos}}, $pos-1, 1;
	# }
	$corrected_peptide = $corrected_peptide . $aa;
	$frame_string = $frame_string . $framer{$pos};
    }
    else {
	my ($previous_frame,$next_frame) = find_prev_next($pos,$q_begin_pos);
	if (defined $previous_frame) {
            $framer{$pos} = $previous_frame;
        }
	elsif (defined $next_frame) {
	    $framer{$pos} = $next_frame;
	}
	else {
	    die "Can't find a previous or upcoming frame:\nPOS:$pos\nBEGIN: $q_begin_pos\nEND:   $q_end_pos\n$read\n\n";
	}
	my $aa = substr $three_frame{$framer{$pos}}, $pos-1, 1;
	$corrected_peptide = $corrected_peptide . $aa;
        $frame_string = $frame_string .$framer{$pos};
        # die "\nERROR: Missing a call at position: $pos\n";
    }
}

print STDOUT ">" . $name . "\n" . $corrected_peptide . "\n";
print STDERR ">" . $name . "\n" . $corrected_peptide . "\n" . $frame_string . "\n";

sub find_prev_next
{
    my $pos = $_[0];
    my $q_begin_pos = $_[1];
    my ($previous_frame,$next_frame);
    for (my $i=1;$i<$pos-$q_begin_pos;$i++) {
	if (exists $framer{$pos-$i}) {
	    $previous_frame = $framer{$pos-$i};
	    last;
	}
    }
    for(my $i=1;$i<$q_end_pos - $pos;$i++) {
	if (exists $framer{$pos+$i}) {
	    $next_frame = $framer{$pos+$i};
	    last;
	}
    }
    return($previous_frame,$next_frame);
}

exit 0;
