package Format;

use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(file_root);
%EXPORT_TAGS = ( DEFAULT => [qw(&fasta_check)]);


sub file_root
{
    my $root = $_[0];
    $root =~ s/.*\///;
    $root = scalar reverse $root;
    $root =~ s/^.*?\.//;
    $root = scalar reverse $root;
    return $root;
}

sub db_root
{
    my $root = $_[0];
    $root =~ s/.*\///;
    return $root;
}

1;
