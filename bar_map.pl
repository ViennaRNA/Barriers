#!/usr/bin/perl -w
# -*-CPerl-*-
# Last changed Time-stamp: <2003-09-02 10:07:48 ivo>

# analyse folding landscape of the growing molecule: 

# read a series of bar files for different length fragments (in the
# order of the growing chain) and compute which minima in successive
# barrier trees are equivalent. Output the correpondence table
# of minima.
# Usage: bar_map.pl  1.bar 2.bar 3.bar ... n.bar

use FindBin qw($Bin);
#use lib ("$Bin", '/home/hofacker/ViennaRNA/Perl/blib/lib', '/home/hofacker/ViennaRNA/Perl/blib/arch');
use RNA;
use RNA::barrier;
use Getopt::Long;
use strict;
use warnings;

Getopt::Long::config("no_ignore_case");

use vars qw/$opt_debug $opt_v $ParamFile $pf $ns_bases/;
&usage() unless GetOptions(
                           "T=f" => \$RNA::temperature,
                           "4" => sub {$RNA::tetra_loop = 0},
                           "d|d0" => sub {$RNA::dangles=0},
                           "d2" => sub {$RNA::dangles=2},,
                           "d3" => sub {$RNA::dangles=3},,
                           "noGU" => \$RNA::noGU,
                           "noCloseGU" => \$RNA::no_closingGU,
                           "noLP" => \$RNA::noLonelyPairs,
                           "logML" => \$RNA::logML,
                           "P=s" => \$ParamFile);

RNA::read_parameter_file($ParamFile) if ($ParamFile);

usage() if $#ARGV<1;
my ($seq1, @lmin1) = read_bar();

my @match_list;
while ($#ARGV >= 0) {
  my ($seq2, @lmin2) = read_bar();

  warn "sequences in bar files don't match"
    unless substr($seq2,0,length($seq1)) eq $seq1;

  my %l2hash = map {$_->[1] => $_->[0]} @lmin2[1..$#lmin2];
  my %match;
  for my $l1 (1..$#lmin1) {
    my $struc = $lmin1[$l1]->[1] . '.' x (length($seq2)-length($seq1));
    my $gs = grad_walk($seq2, $struc);
    if (exists $l2hash{$gs}) {
      $match{$l1} = [$l2hash{$gs}, 0];
    } else {
      my $d=99999999;
      my $best;
      for my $l2 (1..$#lmin2) {
	my $di = RNA::bp_distance($gs, $lmin2[$l2]->[1]);
	($d, $best) = ($di, $l2) if ($di<$d);
      }
      $match{$l1} = [$best, $d];
    }
  }
  $seq1 = $seq2;
  @lmin1 = @lmin2;
  push @match_list, \%match;
}

my @lines = map {[$_]} (1 .. keys(%{$match_list[0]}));
for my $l (0..$#match_list) {
  my %seen;
  my %match  = %{$match_list[$l]};
  for (0..$#lines) {
    my $m = $lines[$_][$l];
    $lines[$_][$l+1] = $match{$m}->[0];
    $seen{$match{$m}->[0]}=1;
#    print "$m -> ", $match{$m}->[0],"\n";
  }
  foreach my $k (sort {$a <=> $b} keys %match) {
    next if $seen{$k};
    my $z = $#lines;
    $lines[$z+1][$l+1] = $k;
  }
#print "\n";
}

@lines = sort {$$a[-1] <=> $$b[-1]} @lines; 
foreach my $l (@lines) {
  print defined($l->[0])? sprintf("%3d", $l->[0]) : "   ";
  for my $b (1..$#{$l}) {
    if (defined($l->[$b])) {
      if ($l->[$b-1] && $match_list[$b-1]{$l->[$b-1]}[1]) {
	print ' ~>';
      } else {
	print ' ->';
      }
      printf("%3d", $l->[$b]);
    }
    else {print "      ";}
  }
  print "\n";
}

sub grad_walk {
  my ($seq, $stru) = @_;
  my $E = RNA::energy_of_struct($seq, $stru);
  my $found;
  do {
    $found = 0;
    foreach my $s (RNA::barrier::get_neighbors($seq, $stru)) {
      my $Ei = RNA::energy_of_struct($seq, $s);
      if ($Ei<$E) {
	($E, $stru) = ($Ei, $s);
	$found++;
      }
    }
  } while($found);
  return $stru;
}


sub read_bar {
  $_ = <>;
  warn "no seq in bar file" unless /^\s+(\S+)/;
  my $seq = $1;
  my @lmin;
  my $nn = 1;
  while (<>) {
    my @F = split;
    next if @F < 5; # shouldn't happen
#    splice(@F,2,1) if ($F[2] eq '('); # got 2 fields from e.g. "( -1.00)"
#    $F[2] =~ s/[()]//g;               # remove brackets from energy field
    $lmin[$nn++] = \@F;
    last if eof;
  }
  return $seq, @lmin;
}

sub usage {
  die "$0 [ViennaRNA options] 1.bar 2.bar [3.bar [...]]\n";
}
