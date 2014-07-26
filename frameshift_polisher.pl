#!/usr/bin/perl -w

# MANUAL FOR frameshift_polisher.pl

=pod

=head1 NAME

frameshift_polisher.pl -- BLAST-driven frameshift remover

=head1 SYNOPSIS

 frameshift_polisher.pl -fasta /Path/to/infile.fasta -outfile /Path/to/output.fasta -db /Path/to/db -work /Path/to/working/directory/
                     [--help] [--manual]

=head1 DESCRIPTION

 BLASTX-drive frameshift polisher. Will perform 6-frame BLAST of sequences against a BLASTable
 database and piece together the frameshifts (if there even are any).
 
=head1 OPTIONS

=over 3

=item B<-f, --fasta>=FILENAME

Input NUCLEOTIDES in FASTA format. (Required) 

=item B<-o, --outfile>=FILENAME

Output FASTA file of PEPTIDEs. (Default: Saved under -work directory as a *.polished.fasta file)

=item B<-d, --db>=FILENAME

Location of the user-defined BLASTable database. (Required) 

=item B<-w, --work>=DIR

Path of working directory, this does not necessarily need to already exist. (Required)

=item B<-t, --threads>=INT

Number of threads to let BLASTX use. (Default = 1)

=item B<-v, --version>

Displays the version. (Optional)

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
use Polisher::Alignment qw(:Both);
use Polisher::Format;

#ARGUMENTS WITH NO DEFAULT
my($fasta,$outfile,$db,$work,$version,$help,$manual);
## Args with defaults
my $threads = 1;

GetOptions (	
				"f|fasta=s"	=>	\$fasta,
                                "o|outfile=s"   =>      \$outfile,
				"d|db=s"	=>	\$db,
                                "w|work=s"      =>      \$work,
                                "t|threads=s"   =>      \$threads,
				"v|version"     =>      \$version,
                                "h|help"	=>	\$help,
				"m|manual"	=>	\$manual);

# VALIDATE ARGS
if ($version) {print " The Frameshift Polisher\n Version 2.0\n Contact: Dan Nasko (dnasko\@udel.edu)\n\n"; exit 0;}
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument -fasta not found.\n\n", -exitval => 2, -verbose => 1)  if (! $fasta );
pod2usage( -msg  => "\n\n ERROR!  Required argument -db not found.\n\n", -exitval => 2, -verbose => 1)  if (! $db );
pod2usage( -msg  => "\n\n ERROR!  Required argument -work not found.\n\n", -exitval => 2, -verbose => 1)  if (! $work );

## QC checks
QC::fasta_check($fasta);
my $num_seqs = QC::nt_check($fasta);
my $infile_root = Format::file_root($fasta);
my $db_root     = Format::db_root($db);
$threads = QC::thread_check($num_seqs, $threads);
unless (defined $outfile) { $outfile = "$work/frameshift_polisher/$infile_root.$db_root.polished.fasta"; }

## Create working directories and update
print `mkdir -p $work`;
print `mkdir -p $work/ncbi-blastx`;
print `mkdir -p $work/frameshift_polisher`;
print " QC checks ... [ Passed ]\n Create working directories ... [ Passed ]\n Begin BLASTx using $threads threads ... ";

## Begin BLASTx
my $blastx_exe = "perl $FindBin::Bin/bin/para_blastx.pl " . 
    "-q $fasta " . 
    "-d $db " . 
    "-o $work/ncbi-blastx/$infile_root.$db_root.btab " . 
    '--outfmt="\"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qseq sseq ppos"\" ' .
    "-t $threads " .
    "-e 1";

print `$blastx_exe`;
print " [BLAST complete]\n";

## Begin parsing the BLASTX output
if ( -z '$work/ncbi-blastx/$infile_root.$db_root.btab') {
    die "\n\n Frameshift Polisher is exiting because none of your sequences found\n a significant hit to any sequences in the BLAST database you provided.\n\n";
}

my $parse_exe = "perl $FindBin::Bin/bin/parse_btab.pl $fasta $work/ncbi-blastx/$infile_root.$db_root.btab" . 
    " > $work/frameshift_polisher/$infile_root.$db_root.rough.fasta" .
    " 2> $work/frameshift_polisher/$infile_root.$db_root.report";
print `$parse_exe`;

my $print_flag = 0;
open(OUT,">$outfile") || die "\n\n Error: Cannot open the output file: $outfile\n\n";
open(IN,"<$work/frameshift_polisher/$infile_root.$db_root.rough.fasta") || die "\n\nError: Unable to open the output fasta from the parser\n\n";
while(<IN>) {
    chomp;
    if ($_ =~ m/^>/) {
	if ($_ =~ m/\[1\]$/) { $print_flag = 1; print OUT "$_\n"; }
    }
    elsif ($print_flag == 1) {
	$print_flag = 0;
	print OUT "$_\n";
    }
}
close(IN);
close(OUT);

print "\n Completed successfully!\n Outputs written to: $outfile\n\n";

exit 0;
