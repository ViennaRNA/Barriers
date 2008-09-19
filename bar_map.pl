#!/usr/bin/perl -w
# -*-Perl-*-
# Last changed Time-stamp: <2006-11-28 15:55:23 xtof>
# $Id: bar_map.pl,v 1.7 2008/09/19 17:17:05 ivo Exp $
#
use FindBin qw($Bin);
use lib "$Bin";
use RNA;
use barrier;
use Getopt::Long;
use Pod::Usage;
use warnings;
use strict;
use vars qw/$opt_debug $opt_v $ParamFile $pf $ns_bases @FILES/;

Getopt::Long::config("no_ignore_case");
pod2usage(-verbose => 0)
  unless GetOptions("T=f"       => \$RNA::temperature,
		    "4"         => sub {$RNA::tetra_loop = 0},
		    "d|d0"      => sub {$RNA::dangles = 0},
		    "d2"        => sub {$RNA::dangles = 2},
		    "d3"        => sub {$RNA::dangles = 3},
		    "noGU"      => \$RNA::noGU,
		    "noCloseGU" => \$RNA::no_closingGU,
		    "noLP"      => \$RNA::noLonelyPairs,
		    "logML"     => \$RNA::logML,
		    "P=s"       => \$ParamFile,
		    "man"       => sub {pod2usage(-verbose => 2)},
		    "help"      => sub {pod2usage(-verbose => 1)});

RNA::read_parameter_file($ParamFile) if ($ParamFile);

pod2usage(-verbose => 0) if $#ARGV<1;
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

# print processed filenames to STDOUT
print "#@FILES\n";

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
  # grow @lines array by adding minima that appear on the next
  # landscape for the first time
  foreach my $k (sort {$a <=> $b} keys %{$match_list[$l+1]}) {
    next if $seen{$k};
    my $z = $#lines;
    $lines[$z+1][$l+1] = $k;
  }
#print "\n";
}

# print correspondence table of local minima to STDOUT
@lines = sort {$$a[-1] <=> $$b[-1]} @lines; 
use Data::Dumper;
#print Dumper(\@lines);
#print Dumper(\@match_list);
foreach my $l (@lines) {
  print defined($l->[0])? sprintf("%3d", $l->[0]) : ' ' x 3;
  for my $b (1..$#{$l}) {
    if (defined($l->[$b])) {
      if ($l->[$b-1] && $match_list[$b-1]{$l->[$b-1]}[1]) {
	print ' ~>';
      } else {
	print ' ->';
      }
      printf("%3d", $l->[$b]);
    }
    else {print ' ' x 6;}
  }
  print "\n";
}

#---
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

#---
sub read_bar {
  $_ = <>;
  print STDERR "Processing $ARGV\n";
  push @FILES, $ARGV;
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

#---
sub usage {
  die "$0 [ViennaRNA options] 1.bar 2.bar [3.bar [...]]\n";
}

=head1 NAME

bar_map.pl - analyse dynamic RNA folding landscapes

=head1 SYNOPSIS

bar_map.pl [-T temp] [-4] [-d[0|1|2|3]] [-no[GU|CloseGU|LP]] [-logML]
           [-P paramfile] 1.bar 2.bar [3.bar [...]]

=head1 DESCRIPTION

The program reads a series of bar files e.g for different length
fragments in the order of the growing chain or for different
temperatures in the order of the heating/cooling schedule, and
computes which minima in successive barrier trees are equivalent. The
program Outputs the correpondence table of minima to I<STDOUT>.

Each column of the output corresponds to the minima of a bar file
(last bar file is the right most gap-less column). The correspondence
between minima of successive columns are calculated by a gradient
walk. Within a row of the output the symbol I<-E<gt>> indicates exact
correspondence between the local minima in the succesive barrier trees
while the symbol I<~E<gt>> indicates approximate correspondence with
maximal similarity (= base pair distance). Approximate similarity
arises if forinstance not all minima are listed in the bar file
because the program C<barriers> was used with the I<-minh> option.

=head1 OPTIONS

=over 4

=item B<-d[0|1|2|3]>

Set the "dangling end" energies for bases adjacent to helices in free
ends and multi-loops: With I<-d1> only unpaired bases can participate
in at most one dangling end. With I<-d2> this check is ignored,
dangling energies will be added for the bases adjacent to a helix on
both sides in any case. With I<-d | -d0> dangling ends are ignored
altogether. With I<-d3> mfe folding will allow coaxial stacking of
adjacent helices in multi-loops. At the moment the implementation will
not allow coaxial stacking of the two interior pairs in a loop of
degree 3 and works only for mfe folding.

=item B<-h, -help>

Display long help message.

=item B<-man>

Display man page.

=item B<-no[GU|CloseGU|LP]>

With I<-noGU> do not allow GU pairs. With I<-noCloseGU> do not allow
GU pairs at the end of helices. With I<-noLP> disallow pairs that can
only occur isolated.

=item B<-P> I<paramfile>

Read energy parameters from  paramfile.

=item B<-T> I<temp>

Rescale energy parameters to a temperature of I<temp>.

=item B<-4>

Turn off special stabilizing energies for certain tetra-loops.

=back

=head1 AUTHORS

Ivo L Hofacker, Peter F Stadler, Christoph Flamm

=head1 BUGS

If in doubt our program is right, nature is at fault. Please send
comments and bug reports to xtof@tbi.univie.ac.at.

=cut

__END__
