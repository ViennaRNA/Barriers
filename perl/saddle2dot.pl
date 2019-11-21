#!/usr/bin/env perl
# -*-Perl-*-
# Last changed Time-stamp: <2004-07-02 13:05:06 xtof>
# $Id: saddle2dot.pl,v 1.1 2004/07/02 11:09:21 xtof Exp $
# Description: convert a barriers-saddel-file in somthing sightly

use Getopt::Long;
use File::Basename;
use strict;
use vars qw/$CUT $NAME %LM/;

#some globals
$CUT = 1000;
$NAME = 'Unknown';
%LM = ();

usage() unless GetOptions("cut=i"  => \$CUT,
			  "name=s" => \$NAME,
			  "h"      => \&usage,
			  "help"   => \&usage);

print_Head();
while (my $line = <>) {
  print_FilledNode($.+1000, 'oval', 'green');
  my @minlist = extract_MinList($line);
  memorize_Mins(\@minlist);
  print_Edge($_, $.+1000) for (@minlist);
}
#print_LabeledNode($_, $_) for (keys %LM);
print_FilledNode($_, 'oval', 'black') for (sort {$a<=>$b} keys %LM);
print_Tail();

#---
sub memorize_Mins { $LM{$_} = 1 for @{shift()} }

#---
sub extract_MinList {
  my $string = shift;
  $string =~ s/^\s+//; # kill leading whitespace
  my @minlist = split(/\s+/, $string);
  splice(@minlist, 0,3); # kill first 3 fields
  return map { $_ <= $CUT ? $_ : () } @minlist;
}

#---
sub print_Edge {
  my($source, $target) = @_;
  return unless @_ == 2;
  print "$source--$target;\n";
}

#---
sub print_FilledNode {
  my($id, $type, $fill) = @_;
  return unless @_ == 3;
  print "$id [shape=circle style=filled fillcolor=$fill];\n";
}

#---
sub print_LabeledNode {
  my($id, $label) = @_;
  return unless @_ == 2;
  print "$id [shape=circle label=$label];\n";
}

#---
sub print_Head { print "graph X {\n" }

#---
sub print_Tail { print "}\n" }

#---
sub usage {
  print STDERR
      "\n  Usage: @{[basename($0)]} [options]\n";
  print STDERR
      "     -cut  <Int>  Set cutoff to <Int>\n",
      "                  (default: $CUT)\n",
      "     -name <Str>  Sets name of graph\n",
      "                  (default: $NAME)\n",
      "     -h           display help message\n\n",
      "  pipe output to a graphviz layout processors (dot, neato, circo)\n",
      "  e.g.   circo -Tps FOO.dot > Foo.ps\n\n";
    exit(0);
}

__END__
