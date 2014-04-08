#!/usr/bin/perl -w

# MANUAL FOR skeleton.pl

=pod

=head1 NAME

frameshift_polisher.pl -- Polishes away frameshifts with BLASTx

=head1 SYNOPSIS

 frameshift_polisher.pl -in /Path/to/infile.fasta -db [RNR] [NR] [POLA]
                     [--help] [--manual]

=head1 DESCRIPTION

 When error correcting PacBio reads it was noticed that many reads still contained a few frameshifts. This pipeline will perform a BLASTx of your nucleotide sequences against either a small marker gene DB or all of NR to find where these frameshifts occur and then correct them.
 
=head1 OPTIONS

=over 3

=item B<-i, --in>=FILENAME

Input file nucleotide sequence file in FASTA format. (Required) 

=item B<-d, --db>=FILENAME

BLASTx database, must choose: RNR or NR (Required)

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
use Getopt::Long;
use File::Basename;
use Pod::Usage;

#ARGUMENTS WITH NO DEFAULT
my($infile,$db,$help,$manual);

GetOptions (	
				"i|in=s"	=>	\$infile,
				"d|db=s"	=>	\$db,
				"h|help"	=>	\$help,
				"m|manual"	=>	\$manual);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument -infile not found.\n\n", -exitval => 2, -verbose => 1)  if (! $infile );
pod2usage( -msg  => "\n\n ERROR!  Required argument -db not found.\n\n", -exitval => 2, -verbose => 1)  if (! $db );

my $database;
if    ($db eq "RNR") { $database = "/home/wommacklab/pipelines/blast_db/uniref90.RNRs.pep"; }
elsif ($db eq "NR")  { $database = "/usr/local/blast_db/nr"; }
elsif ($db eq "POLA"){ $database = "/home/wommacklab/projects/EPIC/marker_genes/pol_A/databases/uniref90_polA"; }
else { die "\n\n Error: The database name you provided is not recognized: $db\n\n Options are 'RNR' or 'NR' or 'POLA'\n\n"; }

print `mkdir -p ~/output_repository`;
print `mkdir -p ~/output_repository/ncbi-blastx`;
print `mkdir -p ~/output_repository/frameshift_polisher`;

my $file_root = $infile;
$file_root =~ s/.*\///;
my $r = scalar reverse $file_root;
$r =~ s/^.*?\.//;
$file_root = scalar reverse $r;
$file_root = $file_root . ".$db";

my $blastx_exe = "blastx " . 
    "-query $infile " . 
    "-db $database " . 
    "-out ~/output_repository/ncbi-blastx/$file_root.btab " . 
    "-outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qseq sseq\" " . 
    "-num_threads 16 " . 
    "-evalue 1e-3";

print `$blastx_exe`;

my $parse_exe = "perl /home/wommacklab/pipelines/bin/frame_shift_corrector.pl " .
    "~/output_repository/ncbi-blastx/$file_root.btab " .
    "> ~/output_repository/frameshift_polisher/$file_root.polished.fasta " .
    "2> ~/output_repository/frameshift_polisher/$file_root.polished.report";

print `$parse_exe`;

my $grep_exe = qq(egrep -A1 " 1\$" ~/output_repository/frameshift_polisher/$file_root.polished.fasta | egrep -v );
$grep_exe = $grep_exe . qw( "\-\-") . qq( >~/output_repository/frameshift_polisher/$file_root.polished.top-hit.fasta);

print `$grep_exe`;


exit 0;
