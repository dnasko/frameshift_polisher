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

QC::fasta_check($fasta);
QC::nt_check($fasta);
my $infile_root = Format::file_root($fasta);
my $db_root     = Format::db_root($db);

print `mkdir -p $work`;
print `mkdir -p $work/ncbi-blastx`;
print `mkdir -p $work/frameshift_polisher`;

print " QC checks ... [ Passed ]\n Create working directories ... [ Passed ]\n Begin BLASTx using $threads threads ... ";

my $blastx_exe = "blastx " .
    "-query $fasta " .
    "-db $db " .
    "-out $work/ncbi-blastx/$infile_root.$db_root.btab" .
    "-outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qseq sseq\" " .
    "-num_threads $threads " .
    "-evalue 1e-3";
print `$blastx_exe`;
print "[complete]\n";



exit 0;
