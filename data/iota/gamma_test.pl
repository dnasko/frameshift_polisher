#!/usr/bin/perl -w
use strict;

my $ref_db     = $ARGV[0];
my $reads_file = $ARGV[1];
my $btab       = $ARGV[2];

my (%Ref,%Reads);
my $header;

open(IN,"<$ref_db") || die "\n cannot open file: $ref_db\n";
while(<IN>) {
    chomp;
    if ($_ =~ m/^>/) {
	$header = $_;
	$header =~ s/^>//;
	$header =~ s/ .*//;
    }
    else {
	$Ref{$header} = $_;
    }
}
close(IN);

open(IN,"<$reads_file") || die "\n cannot open file: $reads_file\n";
while(<IN>) {
    chomp;
    if ($_ =~ m/^>/) {
        $header = $_;
	$header=~ s/^>//;
        $header =~ s/ .*//;
    }
    else {
	$Reads{$header} =$_;
    }
}
close(IN);

my $q_prev = `head -n1 $btab | cut -f1`;
my $r_prev = `head -n1 $btab | cut -f2`;
my $min = 10**9;
my $prev_frame = "";
chomp($q_prev);
open(IN,"<$btab") || die "\n Cannot open $btab\n";
while(<IN>) {
    chomp;
    my @A = split(/\t/, $_);
    if ($A[0] eq $q_prev) {
	if ($min > $A[8]) {
	    $min = $A[8]
	} 
    }
    else {
	my $frame_sign;
	if ($prev_frame > 0) {
	    $frame_sign = "0";
	}
	elsif ($prev_frame < 0) {
	    $frame_sign = "1";
	}
	else {
	    die "Uh... what?\n\n";
	}
	# print STDOUT ">$q_prev\n";
	# print STDERR "perl /Users/dnasko/GitHub/frameshift_polisher/bin/needleman_wunsch__dev.pl -o $min -r $Ref{$r_prev} -q $Reads{$q_prev} -s $frame_sign\n";
	print `perl /Users/dnasko/GitHub/frameshift_polisher/bin/needleman_wunsch__dev.pl -o $min -r $Ref{$r_prev} -q $Reads{$q_prev} -s $frame_sign -n $q_prev`;
	$min= 10**9;
	if ($min > $A[8]) {
            $min = $A[8]
        }
    }
    $q_prev = $A[0];
    $r_prev = $A[1];
    $prev_frame = $A[12];
}
close(IN);

my $frame_sign;
if ($prev_frame > 0) {
    $frame_sign = "0";
}
elsif ($prev_frame < 0) {
    $frame_sign = "1";
}
else {
    die "Uh... what?\n\n";
}
# print ">$q_prev\n";                                                                                                                                               
print `perl /Users/dnasko/GitHub/frameshift_polisher/bin/needleman_wunsch__dev.pl -o $min -r $Ref{$r_prev} -q $Reads{$q_prev} -s $frame_sign -n $q_prev`;

exit 0;
