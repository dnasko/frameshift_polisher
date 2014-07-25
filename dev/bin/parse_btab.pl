#!/usr/bin/perl
use strict;

my $fasta = $ARGV[0];
my $infile = $ARGV[1];

my ($q_prev,$s_prev,$pep);
my $begin = 0;
my $end = 0;
my $line_count = 1;
my $counter = 1;      ## Number of subjects to hit a sequence

my $fixes = 0;        ## Number of frameshift correctiosn for a given query
my $stops = 0;
my $total_fixes = 0;
my $total_stops = 0;

my %Types;

my %Overlaps = (
    "0000" => "Type I",
    "1111" => "Type II",
    "0010" => "Type III",
    "1011" => "Type IV",
    "1010" => "Type V",
    "0011" => "Type VI",
    );

open(IN,"<$infile");
while(<IN>) {
    chomp;
    my @A = split(/\t/, $_);
    my $qid    = $A[0];
    my $sid    = $A[1];
    my $qseq   = $A[14];
    my $sseq   = $A[15];
    my $sstart = $A[8];
    my $send   = $A[9];
    my $num_gaps = $sseq =~ tr/-/-/;
    my $qframe = $A[12];
    my $ppos = $A[16];
    if ($line_count == 1) {
	$fixes += 1;
	$pep = $qseq;
	$begin = $sstart;
	$end = $send;
    }
    else {
	my $binary = "";
	if ($qid eq $q_prev) {
	    if ($sid eq $s_prev && $ppos >= 30) {
		$fixes += 1;
		if ($sstart > $begin) { $binary = "1"; }
		else { $binary = "0"; }
		if ($sstart > $end) { $binary = $binary . "1"; }
		else { $binary = $binary . "0"; }
		if ($send > $begin) { $binary = $binary . "1"; }
                else { $binary = $binary . "0"; }
		if ($send > $end) { $binary = $binary . "1"; }
                else { $binary = $binary . "0"; }
		##############
		## Evaluate ##
		##############
		if ($binary eq "0000") {
		    #print "$line_count: Type I $begin\t$end\t$sstart\t$send\n";
		    $Types{"Type 1"}++;
		    my $diff = $begin - $send;
		    $diff--;
                    for (my $i = 0; $i< $diff; $i++) {
                        $pep = "X" . $pep;
                    }
                    $pep = $qseq . $pep;
                    $begin = $sstart;
		}
		elsif ($binary eq "1111") {
		    #print "$line_count: Type II $begin\t$end\t$sstart\t$send\n";
		    $Types{"Type 2"}++;
		    my $diff = $sstart - $end;
		    $diff--;
		    for (my $i = 0; $i < $diff; $i++) {
			$pep = $pep . "X";
		    }
		    $pep = $pep . $qseq;
		    $end = $send;
		}
		elsif ($binary eq "0010") {
		    #print "$line_count: Type III $begin\t$end\t$sstart\t$send\n";
		    $Types{"Type 3"}++;
		    my $diff = $send - $begin;
		    $diff++;
		    my $to_grab = length($qseq) - $diff;
		    my $sub_seq = substr $qseq, 0, $to_grab;
		    $pep = $sub_seq . $pep;
		    $begin = $sstart;
		}
		elsif ($binary eq "1011") {
		    #print "$line_count: Type IV $begin\t$end\t$sstart\t$send\n";		    
		    $Types{"Type 4"}++;
		    my $diff = $end - $sstart;
                    for (my $i=0; $i <= $diff; $i++) {
			$pep =~ s/.$//;
		    }
		    $pep = $pep . $qseq;
                    $end = $send;
		}
		elsif ($binary eq "1010") {
		    #print "$line_count: Type V $begin\t$end\t$sstart\t$send\n";
		    $Types{"Type 5"}++;
		    my $b_adjust = $sstart - $begin;
		    my $e_adjust = $send - $begin;
		    $b_adjust++;
		    $e_adjust += 2;
		    my $head = substr $pep, 0, $b_adjust;
		    my $tail = substr $pep, $e_adjust;
		    $pep = $head . $qseq . $tail;
		}
		elsif ($binary eq "0011") {
		    #print STDOUT "TYPE6\t$line_count: Type VI $begin\t$end\t$sstart\t$send\n";
		    $Types{"Type 6"}++;
		    my $top_diff  = ($begin - $sstart) + 1;
		    my $tail_diff = ($end - $send) + 1;
		    my $head = substr $qseq, 0, $top_diff;
		    my $tail = substr $qseq, $tail_diff;
		    $pep = $head . $pep . $tail;
		    $begin = $sstart;
		    $end   = $send;
                }
		else {
		    die "\n\n Unexpected overlap class encountered: $binary\n\n";
		}
	    }
	    else {
		if ($counter == 1) {
		    my @A = split(//, $pep);
		    foreach my $i (@A) { if ($i eq "*") {$stops++;}}
		    $pep =~ s/-//g;
		    $pep =~ s/\*/X/g;
		    $fixes--;
		    $total_fixes += $fixes;
		    $total_stops += $stops;
		    print STDERR "$qid [$counter]\t$fixes\t$stops\n";
		    print STDOUT ">$qid" . " [" . $counter . "]\n$pep\n";
		}
		elsif ($counter < 11) {
		    $pep =~ s/-//g;
		    $pep =~ s/\*/X/g;
		    print STDOUT ">$qid" . " [" . $counter . "]\n$pep\n";
		}
		$fixes = 0;
		$stops = 0;
		$pep = $qseq;
		$begin = $sstart;
		$end = $send;
		$counter++;
	    }
	}
	else {
	    if ($counter < 11) {
		$pep =~ s/-//g;
		$pep =~ s/\*/X/g;
		print STDOUT ">$q_prev" . " [" . $counter . "]\n$pep\n";
	    }
	    $pep = $qseq;
	    $begin = $sstart;
	    $end = $send;
	    $counter = 1;
	    $fixes++;
	}
    }
    $q_prev = $qid;
    $s_prev= $sid;
    $line_count++;
}
close(IN);

print STDERR "\nTOTAL\tFrameshifts\tStop Codons
\t$total_fixes\t\t$total_stops\n\n";

foreach my $i (sort keys %Types) {
    print STDERR "$i\t$Types{$i}\n";
}

exit 0;
