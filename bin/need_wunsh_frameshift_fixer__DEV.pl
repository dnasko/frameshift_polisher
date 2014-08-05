#!/usr/bin/perl -w

# MANUAL FOR skeleton.pl

=pod

=head1 NAME

need_wunsh_frameshift_fixer.pl -- frameshift correction with iterative Needleman-Wunsch alignments

=head1 SYNOPSIS

 need_wunsh_frameshift_fixer.pl -query ATGCTAG -ref MLCKWTGG -offset 245 [-antisense]
                     [--help] [--manual]

=head1 DESCRIPTION

 Blah.
 
=head1 OPTIONS

=over 3

=item B<-q, --query>=ATGCTA

Nucleotide query sequence. (Required) 

=item B<-r, --ref>=MLKCCTW

Peptide reference sequence. (Required) 

=item B<-o, --offset>=INT

Reference offset. (Required)

=item B<-a, --antisense>

Perform a 3 frame translation in the antisense. (Optional)

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
use FindBin;
use Cwd 'abs_path';
use lib abs_path("$FindBin::Bin/../lib");
use File::Basename;
use Polisher::QC qw(:Both);
use Polisher::Alignment qw(:Both);
use Polisher::Format;

#ARGUMENTS WITH NO DEFAULT
my($query,$ref,$offset,$antisense,$help,$manual);

GetOptions (	
			 "q|query=s"	=>	\$query,
                         "r|ref=s"      =>      \$ref,
                         "o|offset=i"	=>	\$offset,
                         "a|antisense"  =>      \$antisense,
			 "h|help"	=>	\$help,
			 "m|manual"	=>	\$manual);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument -query not found.\n\n", -exitval => 2, -verbose => 1)  if (! $query );
pod2usage( -msg  => "\n\n ERROR!  Required argument -ref not found.\n\n", -exitval => 2, -verbose => 1)  if (! $ref );
pod2usage( -msg  => "\n\n ERROR!  Required argument -offset not found.\n\n", -exitval => 2, -verbose => 1)  if (! $offset );

my ($frame_start,$frame_end) = (1,3);
if ($antisense) { ($frame_start,$frame_end) = (4,6); }

for (my $frame=$frame_start; $frame <= $frame_end; $frame++) {
    my $query_peptide = Alignment::translate($query, $frame);
    print "$frame\n$query_peptide\n\n";
}

exit 0;
