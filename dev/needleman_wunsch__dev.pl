#!/usr/bin/perl -w

use strict;
use FindBin;
use Cwd 'abs_path';
use lib abs_path("$FindBin::Bin/lib");
use File::Basename;
use Polisher::QC qw(:Both);
use Polisher::Alignment qw(:Both);
use Polisher::Format;

my $reference_offset = $ARGV[0];
my $reference_file   = $ARGV[1];
my $read_file        = $ARGV[2];

my ($reference,$read);

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

my $trmd_reference = substr $reference, ($reference_offset-1);
my %three_frame;

for (my $frame=1; $frame <=3; $frame++) {
    my $query_peptide = Alignment::translate($read, $frame);
    $three_frame{$frame} = $query_peptide;
    print "Frame: $frame\n$query_peptide\n\n";
}


exit 0;
