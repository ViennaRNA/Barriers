#!/usr/bin/perl
# -*-Perl-*-
# Last changed Time-stamp: <2006-07-19 12:29:41 xtof>
# $Id: to_bar_in.pl,v 1.1 2006/07/25 14:34:32 xtof Exp $

use strict;

my @DATA = ();

while (<>) {
  my @F = split;
  push @DATA, [@F];
}

@DATA = sort {   $a->[1]<=>$b->[1]
	      || $a->[2]<=>$b->[2]
	      || $a->[3]<=>$b->[3]
	      || $a->[4]<=>$b->[4] } @DATA;

print join(" ", @$_), "\n" for @DATA;


__END__
