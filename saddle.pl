#!/usr/bin/perl -w
# -*-Perl-*-
# Last changed Time-stamp: <2001-01-25 10:07:25 ivo>

use RNA;
use Getopt::Long;
use strict;

 Getopt::Long::config('no_ignore_case');

use vars '$max', '$max_en', '$opt_4', '$opt_d', '$opt_d2', '$ParamFile',
    '$opt_v','$opt_d3' ;
$max = 20;
$max_en = 99999.;

&usage() unless GetOptions("T=f" => \$RNA::temperature,
                           "4", "d", "d2", "d3",
                           "logML" => \$RNA::logML,
                           "P=s" => \$ParamFile,
                           "max=i"=> \$max,
			   "E=f" => \$max_en,
                           "v");

$RNA::tetra_loop = 0 if ($opt_4);
$RNA::dangles = 0 if ($opt_d);
$RNA::dangles = 2 if ($opt_d2);
$RNA::dangles = 3 if ($opt_d3);

 RNA::read_parameter_file($ParamFile) if (defined($ParamFile));

# importing barrier.pl will call RNA::update_fold_params()
# this should really be  use RNA::barrier;
use lib '/home/linse/ivo/RNA/barriers';
use barrier;
 
 RNA::read_parameter_file($ParamFile) if (defined($ParamFile));
RNA::update_fold_params();
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
     RNA::barrier::find_saddle($string, $struc1, $struc2, $max, $max_en);
   $saddle = "not found" if (!defined($saddle));
   if (($opt_v)&&($#RNA::barrier::path>-1)) {
      foreach my $s (@RNA::barrier::path) {
	 printf "$s (%6.2f)\n", RNA::energy_of_struct($string, $s);
      }
      print "\n";
   } else {
      printf "$saddle (%6.2f)\n", $saddleE;
   }
}

sub usage {
   die "$0 [-T] [-4] [-d[2]] [-logML] [-P file] [-max m] [-E maxE]\n",
   "-max m \tlimit number of paths at each depth to at most m\n",
   "-E max \tsearch for saddle points with energy < maxE\n" .
       "\tstandard options for RNA energy parameters\n";
}
# End of file
