#!/usr/bin/perl
# -*-Perl-*-
# Last changed Time-stamp: <2006-07-18 15:35:54 xtof>
# $Id: to_bar_in.pl,v 1.1 2006/07/18 14:10:07 xtof Exp $

use Getopt::Long;
use POSIX;
use strict;
use vars qw/@DATA $COLS/;

my $COLS = 4;

usage() unless GetOptions("cols=i" => \$COLS);

my $params = read_data();
my $stdv   = calculate_stdv(\@DATA);

# normalize data
for my $line (@DATA) {
  my $sum = 0;
  for (2..$COLS-1) {
#    $line->[$_] /= $stdv->[$_];
#    $sum += $line->[$_]/= $stdv->[$_];
#    $line->[$_] *= 1000 if $_ == 2;
    $line->[$_] = POSIX::ceil($line->[$_]);
    $sum += ($line->[$_]/ $stdv->[$_]);
  }
  push @$line, $sum;
}

# poset sort data
@DATA = sort {   $a->[1] <=> $b->[1]
	      || $b->[2] <=> $a->[2]
              || $b->[3] <=> $a->[3]
	      || $b->[4] <=> $a->[4]} @DATA;

# output barriers input file
print to_header(@$params, $COLS-1), "\n";
printf "%s %6.2f %4d %4d %4d\n", @$_[0..$COLS] for @DATA;

#---
sub read_data {
  my ($structure, $constrain, $protein, $max_m, $min_as);
  while (<>) {
    # ignore lines
    next if m/^\*/;
    next if m/^\s+/;

    # memorize parameters
    $structure = trim_whitespace((split /:/, $_, 2)[1]), next if m/^Structure/;
    $constrain = trim_whitespace((split /:/, $_, 2)[1]), next if m/^Constrain/;
    $protein   = trim_whitespace((split /:/, $_, 2)[1]), next if m/^Protein/;
    $max_m     = trim_whitespace((split /:/, $_, 2)[1]), next if m/^Max/;
    $min_as    = trim_whitespace((split /:/, $_, 2)[1]), next if m/^Min/;
    
    # ignore lines
    next if m/:/;
    
    # memorize data
    push @DATA, [(split)[0..$COLS-1]];
#    $DATA[-1]->[2] *= 10000; # scale frequency in ensemble
    $DATA[-1]->[2] = exp -30 if $DATA[-1]->[2] == 0;
    $DATA[-1]->[2] = log $DATA[-1]->[2];
  }

  return [$structure, $constrain, $protein, $max_m, $min_as]
}

#---
sub calculate_stdv {
  my (@S, @S2, @stdv);

  my $data = shift;
  my $N = $#$data+1;
  my $n = $#{$data->[0]};
  my @S  = map { 0 } 0 ..$n;
  my @S2 = map { 0 } 0 ..$n;

  for my $line (@$data) {
    for (my $i = 0; $i <= $n; $i++) {
      $S[$i]  += $line->[$i];
      $S2[$i] += $line->[$i] * $line->[$i];
    }
  }

  for (my $i=0; $i <= $n; $i++) {
    push @stdv, sqrt(($S2[$i] - $S[$i]*$S[$i]/$N)/($N));
  }

  return \@stdv;
}

#---
sub trim_whitespace {
  local $_ = shift;
  s/\s+//g;

  return $_;
}

#---
sub to_header {
  my ($structure, $constrain, $protein, $max_m, $min_as, $cols) = @_;
  my $header =
      join " ", $constrain, -999999, ":$cols:", "SECIS,$max_m,$min_as";

  return join "\n", $header, join " ", $structure, $protein;
}

#---
sub usage {
  printf STDERR "\nusage: @{[basename $0]} [--cols number]\n";
  printf STDERR "default: cols = $COLS\n";
  exit;
}


__END__
