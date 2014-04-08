#!/usr/bin/perl -w

# MANUAL FOR frameshift_polisher.pl

=pod

=head1 NAME

frameshift_polisher.pl -- BLAST-driven frameshift remover

=head1 SYNOPSIS

 frameshift_polisher.pl -fasta /Path/to/infile.fasta -db /Path/to/db -work /Path/to/working/directory/
                     [--help] [--manual]

=head1 DESCRIPTION

 BLASTX-drive frameshift polisher. Will perform 6-frame BLAST of sequences against a BLASTable
 database and piece together the frameshifts (if there even are any).
 
=head1 OPTIONS

=over 3

=item B<-f, --fasta>=FILENAME

Input NUCLEOTIDES in FASTA format. (Required) 

=item B<-d, --db>=FILENAME

Location of the user-defined BLASTable database. (Required) 

=item B<-w, --work>=DIR

Path of working directory, this does not necessarily need to already exist. (Required)

=item B<-t, --threads>=INT

Number of threads to let BLASTX use. (Default = 1)

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-m, --manual>

Displays full manual.  (Optional) 

=back

=head1 DEPENDENCIES

Requires the following Perl libraries.



=head1 AUTHOR

Written by Daniel Nasko, 
Center for Bioinformatics and Computational Biology, University of Delaware.

=head1 REPORTING BUGS

Report bugs to dnasko@udel.edu

=head1 COPYRIGHT

Copyright 2014 Daniel Nasko.  
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
This is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  

Please acknowledge author and affiliation in published work arising from this script's 
usage <http://bioinformatics.udel.edu/Core/Acknowledge>.

=cut


use strict;
use FindBin;
use Cwd 'abs_path';
use lib abs_path("$FindBin::Bin/lib");
use Getopt::Long;
use File::Basename;
use Pod::Usage;
use Polisher::QC qw(:Both);
use Polisher::Format;

#ARGUMENTS WITH NO DEFAULT
my($fasta,$db,$work,$help,$manual);
## Args with defaults
my $threads = 1;

GetOptions (	
				"f|fasta=s"	=>	\$fasta,
				"d|db=s"	=>	\$db,
                                "w|work=s"      =>      \$work,
                                "t|threads=s"   =>      \$threads,
				"h|help"	=>	\$help,
				"m|manual"	=>	\$manual);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument -fasta not found.\n\n", -exitval => 2, -verbose => 1)  if (! $fasta );
pod2usage( -msg  => "\n\n ERROR!  Required argument -db not found.\n\n", -exitval => 2, -verbose => 1)  if (! $db );
pod2usage( -msg  => "\n\n ERROR!  Required argument -work not found.\n\n", -exitval => 2, -verbose => 1)  if (! $work );

## QC checks
QC::fasta_check($fasta);
QC::nt_check($fasta);
my $infile_root = Format::file_root($fasta);
my $db_root     = Format::db_root($db);

## Create working directories and update
print `mkdir -p $work`;
print `mkdir -p $work/ncbi-blastx`;
print `mkdir -p $work/frameshift_polisher`;
print " QC checks ... [ Passed ]\n Create working directories ... [ Passed ]\n Begin BLASTx using $threads threads ... ";

## Begin BLASTX
my $blastx_exe = "blastx " .
    "-query $fasta " .
    "-db $db " .
    "-out $work/ncbi-blastx/$infile_root.$db_root.btab" .
    "-outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qseq sseq\" " .
    "-num_threads $threads " .
    "-evalue 1e-3";
print `$blastx_exe`;
print "[complete]\n";

## Begin parsing the BLASTX output
if ( -z '$work/ncbi-blastx/$infile_root.$db_root.btab') { die "\n\n Frameshift Polisher is exiting because none of your sequences found a significant hit to any sequences in the BLAST database you provided.\n\n" }

my ($q_prev,$s_prev,$pep);

my $line_count = 1;
my $begin = 0;
my $end = 0;
my $counter = 1;
my %Types;

my $fixes = 0;
my $stops = 0;
my $total_fixes = 0;
my $total_stops = 0;

my %Overlaps = (
    "0000" => "Type I",
    "1111" => "Type II",
    "0010" => "Type III",
    "1011" => "Type IV",
    "1010" => "Type V",
    "0011" => "Type VI",
    );

open(IN,"<$work/ncbi-blastx/$infile_root.$db_root.btab") || die "\n\n Error: Cannot open BLAST's output: $work/ncbi-blastx/$infile_root.$db_root.btab\n\n";
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
    if ($line_count == 1) {
	$fixes += 1;
	$pep = $qseq;
	$begin = $sstart;
	$end = $send;
	#print "$begin\t$end\t$pep\n";
    }
    else {
	my $binary = "";
	if ($qid eq $q_prev) {
	    if ($sid eq $s_prev) {
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
                    for (my $i = 0; $i < $diff; $i++) {
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
		    my $diff = $sstart - $send;
                    $diff++;
                    my $sub_seq = substr $qseq, $diff;
                    $pep = $pep . $sub_seq;
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
		        #print "$line_count: Type VI $begin\t$end\t$sstart\t$send\n";
		    $Types{"Type 6"}++;
		        
                }
		else {
		    die "\n\n Unexpected overlap class encountered: $binary\n\n";
		}
	    }
	    else {
		if ($counter == 1) {
		    my @A = split(//, $_);
		    foreach my $i (@A) { if ($i eq "*") {$stops++;}}
		    $pep =~ s/-//g;
		    $pep =~ s/\*/X/g;
		    $fixes--;
		    $total_fixes += $fixes;
		    $total_stops += $stops;
		    print STDERR "$qid -> $counter\t$fixes\t$stops\n";
		}
		print STDOUT ">$qid" . " -> " . $counter . "\n$pep\n";
		$fixes = 0;
		$stops = 0;
		$pep = $qseq;
		$begin = $sstart;
		$end = $send;
		$counter++;
	    }
	}
	else {
	    $pep =~ s/-//g;
	    $pep =~ s/\*/X/g;
	    print STDOUT ">$q_prev" . " -> " . $counter . "\n$pep\n";
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
