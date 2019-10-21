#!/usr/bin/perl -w
# -*-Perl-*-
# Last changed Time-stamp: <2005-10-14 14:33:53 ivo>
# $Id: treeplot.pl,v 1.2 2005/10/14 12:38:27 ivo Exp $

use Getopt::Long;
use File::Basename;
use strict;
use vars qw/$DATA @BBOX $TITLE $L2I $I2L $A $B/;

use constant {
  ID       => 0, # uniq label of celltype
  ENERGY   => 2, # deathtime  of celltype
  FATHER   => 3, # uniq label of father celltype
  LIFETIME => 4, # lifetime of celltype (birthtime - deathtime)
  SADDEL   => 5, # saddle point hight
  IDX      => 6, # line index in input file
};

# some globals
$TITLE = $ARGV[0] || "stdin"; # title of ps-file
@BBOX  = (72, 144, 522, 700); # bounding box
$DATA  = {}; # hash holding input data
$L2I   = {}; # label => index  hash
$I2L   = {}; # index => label  hash

parse_phylo_file();
idx_list_sorted_by_SADDEL();
x_order_leafs();
print_eps();

#---
# sub parse_bar_file {
#   local $_;
#   $_ = <>;
#   while (<>) {
#     my @F = split [0,2,3,4]; # ID struct Energy Father barrier
#     $DATA->{$F[ID]} = [@F, $F[ENERGY]+$F[LIFETIME], $.-1];
#     $L2I->{$F[ID]} = $.-1;
#     $I2L->{$.-1}   = $F[ID];
#   }
# }
sub parse_phylo_file {
  local $_;
  $_ = <>;
  while (<>) {
    my @F = split; # ID Energy Father LifeSpan
#    $F[4] += 0.001 if $F[3]==0;
    $DATA->{$F[ID]} = [@F, $F[ENERGY]+$F[LIFETIME], $.-2];
    $L2I->{$F[ID]} = $.-2;
    $I2L->{$.-2}   = $F[ID];
  }
#  $L2I->{0} = 0;
#  $I2L->{0} = 0;
}

#---
sub x_order_leafs {
  use constant { VAL => 0, NXT => 1 };
  my ($l, $i, $f);
  # an ugly mimicry of a linked list jak!
  my @chain = map {[$_,undef,undef]} $[..$#$A;
  for my $k (map{$A->[$_]} $[..$#$A) {
    $f = $L2I->{$DATA->{$I2L->{$k}}->[FATHER]};
    $f = 0 unless defined $f;
#    print STDERR "$k $f\t";
    next if $k == $f; # lowest node doesn't merge
    for ($l=$chain[$f]; defined($l->[NXT]); $l=$chain[$l->[NXT]]) {}
    $l->[NXT] = $k; # attach child to chain of father
  }

  # chain[$f] now starts the orderd chain, next fill in the num field
  for ($i=0, $l = $chain[$f];
       defined($l->[NXT]); $l=$chain[$l->[NXT]]) { $l->[VAL] = $i++ }
  $l->[VAL] = $i;
  $B = [map{$chain[$_]->[VAL]} $[..$#$A];
}

#---
sub print_eps {
  print_eps_head();
  print_label_array();
  print_leaf_array();
  print_saddle_array();
  print_eps_body_and_tail();
}

#---
sub print_label_array {
  my $i = 0;
  print "  /LABEL [";
  for my $label (map{$I2L->{$_}} $[..$#$A) {
    ($i++ % 10) ? print " " : print "\n   ";
    printf "%g", $label;
  }
  printf "  ] def\n";
}

#---
sub print_leaf_array {
  my $i = 0;
  printf "%% leaf node coordinates\n";
  printf "  /LEAF [";
  for my $idx ($[..$#$B) {
    ($i++ % 5) ? print " " : print "\n   ";
    printf
	"[%-3d %7.3f]",
	$B->[$idx],
	$DATA->{$I2L->{$idx}}->[ENERGY];
  }
  printf "  ] def\n";
}

#---
sub print_saddle_array {
  my $i = 0;
  printf "%% internal nodes (saddle) coordinates, sorted by height\n";
  printf "  /SADDEL [";
  for my $idx (@$A) {
    ($i++ % 4) ? print " " : print "\n   ";
    my $f = $L2I->{$DATA->{$I2L->{$idx}}->[FATHER]};
    $f = -1 if !defined($f);
    printf
	"[%3d %3d %7.3f]",
	$idx,
	$f,
	$DATA->{$I2L->{$idx}}->[SADDEL];
  }
  printf "  ] def\n";
}

#---
sub idx_list_sorted_by_SADDEL {
  $A = [ map{$DATA->{$_}->[IDX]} (sort by_SADDEL keys %$DATA) ];
}

#---
sub by_SADDEL { int(1000*$DATA->{$a}->[SADDEL]+0.1) <=> 
		int(1000*$DATA->{$b}->[SADDEL]+0.1) ||
		$DATA->{$b}->[FATHER] <=> $DATA->{$a}->[FATHER] ||
		$DATA->{$a}->[ENERGY] <=> $DATA->{$b}->[ENERGY]}

#---
sub print_eps_head {
  print << "EOH";
%!PS-Adobe-2.0 EPSF-1.2
%%Title: $TITLE
%%Creator: @{[basename($0)]}
%%CreationDate: @{[scalar(localtime())]}
%%Pages: 1
%%BoundingBox: $BBOX[0] $BBOX[1] $BBOX[2] $BBOX[3]
%%EndComments
%%BeginProlog
/treedict 100 dict def
treedict begin
% x y  => min(x,y)
  /min { 2 copy gt { exch } if pop } bind def
  /max { 2 copy lt { exch } if pop } bind def
  /cmtx matrix currentmatrix def
  /STR 128 string def
  /NumH 1 def
% - => -
  /Init {
    /LX [
      LEAF {0 get} forall
    ] def

    /Helvetica findfont fsize scalefont setfont
    /Lo [
      (X) stringwidth pop % width
      newpath 0 0 moveto
      (X) true charpath
      flattenpath pathbbox
      pop exch pop exch sub neg 2 div % height
     ] def
  } def
% - => -
  /DrawScale {
  gsave 
    maxy miny sub 30 div dup maxy add /maxy exch def miny sub /miny def
    maxy miny sub log 0.9 sub floor 10 exch exp /tick exch def
    newpath
    LEAF length 0.5 sub 0 translate 0 miny moveto 0 maxy miny sub rlineto
    miny tick div ceiling tick mul dup 0 exch moveto 
    maxy exch sub tick div cvi 1 add dup { % draw minor ticks
      0.15 0 rlineto
      -0.15 tick rmoveto
    } repeat
    % calculate major tick spacing (10, 5, or 2 minor ticks)
    dup 69 gt { pop 10
    } {
      32 gt { 5 }
      {2} ifelse
    } ifelse
    tick mul /mtick exch def
    miny mtick div ceiling mtick mul dup 0 exch moveto
    maxy exch sub mtick div cvi 1 add {
      0.3 0 rlineto 
      gsave currentpoint 10 mul round 10 div cmtx setmatrix
      STR cvs dup stringwidth pop
      Lo aload pop 3 1 roll add neg exch rmoveto show pop
      grestore
      -0.3 mtick rmoveto
    } repeat
    cmtx setmatrix stroke    
  grestore
  } def
% - => -
  /SetBarFont {
    matrix currentmatrix cmtx setmatrix
    /Helvetica findfont fbsize scalefont setfont
    setmatrix
  } bind def
% - => -
  /SetLabelFont {
    matrix currentmatrix cmtx setmatrix
    /Helvetica findfont fsize scalefont setfont
    setmatrix
  } bind def
% str => -
  /Rotshow {
    gsave
      cmtx setmatrix -90 rotate
      Lo aload pop
      rmoveto show
    grestore
  } def
% dy => - 
  /Rlineto {
    dup abs MinHeight ge { % draw height at middle of line
      dup gsave
	dup 2 div 0 exch rmoveto
	cmtx setmatrix -90 rotate
	abs STR cvs dup stringwidth pop 2 div neg
	//NumH rmoveto
	show
      grestore
    } if
    0 exch rlineto
  } def
% - => -
  /Drawlabels {
   0 LEAF {
      aload pop moveto
      dup LABEL exch get STR cvs Rotshow
      1 add
    } forall pop
  } def
% n => n\'    Detect whether a minimum is connected
  /MRX {
     /murxi { true } def
     dup 0 lt { pop 0 /murxi { false } def } if
  } def
% - => -
  /Connectlmins {
    newpath
    SADDEL {
      /forest {false} def  %  draw as tree or forest node
      aload pop exch dup 0 lt { pop 0 /forest {true} def} if   % => c h f
      dup LX exch get [ exch LX 5 index get add 2 div % => c h f [ nx
      3 index ]				         % => c h f [ nx h ]
      3 -1 roll dup LEAF 6 -1 roll get aload pop % => f [nx h] h h cx cy
      dup 3 1 roll moveto		         % => f [] h h cy
      sub Rlineto                                % => f [] h
      LEAF 3 index get aload pop exch		 % => f [] h fy fx
      2 index forest {moveto} {lineto} ifelse 
      sub neg Rlineto			         % => f [] h fy
      LEAF 3 1 roll put
    } forall
    gsave
      cmtx setmatrix stroke
    grestore
  } def
end
%%EndProlog
treedict begin
% data starts here!!!
EOH
}

#---
sub print_eps_body_and_tail {
  print << "EOB";
  /fsize 10 def
  /fbsize 7 def
  /yoffset 1 def
  Init
  521 144 yoffset mul fsize 1.5 mul add translate
%  1.0 0.8 scale
  72 521 sub LEAF length div % x-scale
  699 144 fsize dup add add sub
  SADDEL dup length 1 sub get 2 get /maxy exch def % max height
  9999999 LEAF { aload pop exch pop min } forall
  /miny exch def % min height
  maxy miny sub dup 20 div /MinHeight exch def
%  /MinHeight { 99999999 } def	
  div scale
  .5 miny neg translate
  SetLabelFont
  Drawlabels
  DrawScale
  SetBarFont
  Connectlmins
  showpage
end
%%EOF  
EOB
}

#---
sub usage {
  print STDERR
      "\n  Usage: @{[basename($0)]} [options] FOO\n",
      "     -h           display help message\n\n",
      "  FOO must numerically be sorted by 2nd column (sort +1n FOO)\n",
      "  Line-Format of FOO must be:\n",
      "  Id -DeathTime FathersId DeathTime-Birthtime\n\n";
  exit(0);
}



__END__
