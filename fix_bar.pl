#!/usr/local/bin/perl
# -*-Perl-*-
# Last changed Time-stamp: <1999-07-30 00:08:33 ivo>

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


# read in old .bar file
my $first_line = <>;
my $sequence = $1 if ($first_line =~ /(\S+)/);
my @lmin;
my %unknown;
my $print_saddle=0;

while (<>) {
    my @F = split;
    splice(@F,2,1) if ($F[2] eq '('); # got 2 fields from e.g. "( -1.00)"
    $F[2] =~ s/[()]//g;               # remove brackets from energy field
    $F[4] += $F[2];         # use saddle energy instead of barrier height
    $print_saddle = 1 if $F[5] =~ /\.\./;  # file contains saddle structures
    my $nn = shift @F;
    $lmin[$nn] = \@F;
    $unknown{$nn} = '' if (($F[2]==0)&&($nn>1));  # father = 0
}

foreach my $nn (sort {$a <=> $b} keys %unknown) {
    my ($struc, $en, $father, $sE, @rest) = @{$lmin[$nn]};
    
    $sE = 9999.;  # will contain saddle energy
    my $ss;          # saddle structure
    
    foreach my $l (1..$nn-1) {    # compare with other lmins of lower energy
	my ($struc2, $en2, $f, @rest) = @{$lmin[$l]};

	my ($saddleE, $saddle) = 
	  RNA::barrier::find_saddle($sequence, $struc, $struc2, $max, $sE);
	if (defined($saddle)) { # found something
	    my $ff = $l;
	    while (($f!=0)&&($saddleE>$lmin[$ff]->[3])) {
                # find the father's father ....
		$ff = $f;
		$f = $lmin[$f]->[2];
		if (($sE<($lmin[$ff]->[3]))&&
		    ($father<$ff)) { # $ff has incorrect saddle!
		    print STDERR "$nn $l $sE $lmin[$ff]->[3] $f $ff $saddleE\n";
		    warn "inconsitent energies" if (!exists($unknown{$ff}));
#		    print "$ff $lmin[$ff]->[3] $father $sE\n";
		    $lmin[$ff]->[3] = $sE;
		    $lmin[$ff]->[2] = $father;
		    $lmin[$ff]->[4] = $ss if ($print_saddle);
		}
	    }
	    ($sE, $ss, $father) = ($saddleE, $saddle, $ff);
	}
    }
    # store lowest barrier
    if ($print_saddle) {
	$lmin[$nn] = [$struc, $en, $father, $sE, $ss, @rest];
    } else {
	$lmin[$nn] = [$struc, $en, $father, $sE, @rest];
    }
}


# write the new bar file output
print $first_line;
foreach my $nn (1..$#lmin) {
    $lmin[$nn]->[3] -=  $lmin[$nn]->[1];
    printf "%4d %s (%6.2f) %4d %6.2f %s\n", $nn, @{$lmin[$nn]};
}
    

sub usage {
    die "$0 [-T] [-4] [-d[2]] [-logML] [-P file] [-max m] file.bar\n",
    "-max m \tlimit number of paths at each depth to at most m\n",
    "\tstandard options for RNA energy parameters\n",
    "adds missing barrier heights (approximate) to file.bar\n";   
}

# End of file
