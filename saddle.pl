#!/usr/local/bin/perl -w
# -*-Perl-*-
# Last changed Time-stamp: <1999-07-29 20:40:44 ivo>

use RNA;
use Getopt::Long;
use strict;

 Getopt::Long::config('no_ignore_case');

use vars '$max', '$opt_4', '$opt_d', '$opt_d2', '$ParamFile', '$opt_v';
$max = 20;

&usage() unless GetOptions("T=f" => \$RNA::temperature,
                           "4", "d", "d2",
                           "logML" => \$RNA::logML,
                           "P=s" => \$ParamFile,
                           "max=i"=> \$max,
                           "v");

$RNA::tetra_loop = 0 if ($opt_4);
$RNA::dangles = 0 if ($opt_d);
$RNA::dangles = 2 if ($opt_d2);

require 'barrier.pl';
    
print "Input sequence and two structures\n" .
    "....,....1....,....2....,....3....,....4" .
    "....,....5....,....6....,....7....,....8\n"
    if ((-t STDIN) && (-t STDOUT) && ($#ARGV<0));

$, = ' ';

while (<>) {
    my $string;
    
    next if (/^>\s*(\S+)/);

    last if (/^@/);
    if (/(\S+)/) {
        $string = $1;
    } else {
        next;
    }
    $string = uc($string);
    my $length = length($string);
    $_ = <>;
    my $struc1 = $1 if (/([(.)]+)/);
    die ("unequal length: $string $struc1\n$_\n")
        unless (length($struc1)==$length);
    $_ = <>;
    my $struc2 = $1 if (/([(.)]+)/);
    die ("unequal length: $string $struc2\n$_\n")
	unless (length($struc2)==$length);

    my ($saddleE, $saddle) =
      RNA::barrier::find_saddle($string, $struc1, $struc2, $max);

    printf "$saddle (%6.2f)\n", $saddleE;
}

sub usage {
    die "$0 [-T] [-4] [-d[2]] [-logML] [-P file] [-max m]\n" .
	"-max m \tlimit number of paths at each depth to at most m\n" .
	    "\tstandard options for RNA energy parameters\n";
}
# End of file
