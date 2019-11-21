#!/usr/bin/env perl
# -*-Perl-*-
# Last changed Time-stamp: <2001-07-11 21:58:05 xtof>

use Getopt::Long;
use strict;
use vars qw/$CUT %LM/;
$CUT = 1000;
%LM = ();
usage() unless GetOptions("cut=i" => \$CUT);
print STDERR "MESSAGE: cut is $CUT\n";
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
sub memorize_Mins {
  my $minlist = shift;
  $LM{$_} = 1 for @$minlist;
}

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
  print
      "edge [\n",
      "source $source\n",
      "target $target\n",
      "]\n";
	  
#---
sub print_FilledNode {
  my($id, $type, $fill) = @_;
  return unless @_ == 3;
  print
      "node [\n",
      "id $id\n",
      "graphics [\n",
      "type \"$type\"\n",
      "fill \"$fill\"\n",
      "]\n",
      "]\n";
}

#---
sub print_LabeledNode {
  my($id, $label) = @_;
  return unless @_ == 2;
  print
      "node [\n",
      "id $id\n",
      "label \" $label\"\n",
      "]\n";
}

#---
sub print_Head {
  print
      "graph [\n",
      "version 2\n",
      "directed 0\n",
      "node_style [\n",
      "name \"default_node_style\"\n",
      "style [\n",
      "graphics [\n",
      "w 16.0\n",
      "h 16.0\n",
      "]\n",
      "]\n",
      "]\n",
      "edge_style [\n",
      "name \"default_edge_style\"\n",
      "style [\n",
      "graphics [\n",
      "]\n",
      "]\n",
      "]\n";
}

#---
sub print_Tail { print "]\n" }
      
	  
  
}

__END__
