#!/usr/bin/perl
# -*-CPerl-*-
# Last changed Time-stamp: <2017-10-12 18:42:32 mtw>

use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use File::Basename;
use strict;
use warnings;

my ($sequence,$unbound,$bound,$fhu,$fhb);
my ($basename,$bardir,$suffix);
my $bar=undef;
my $rates=undef;
my $binrates=undef;
my %lmins_b = (); # indices of ligand-bound lmins
my %lmins_u = (); # indices of unbound lmins
my @bar_o   = (); # AoH holding the original bar file

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions(
					   "b|bar=s"    => \$bar,
					   "l|log=s"    => \&set_logfile,
					   "man"        => sub{pod2usage(-verbose => 2)},
                                           "h|help"     => sub{pod2usage(1)}
					  );

unless (-f $bar){
  warn "Could not find input BARRIERS input file provided via -b|--bar option";
  pod2usage(-verbose => 0);
}
if($bar =~ m/\.bar$/){
  ($basename,$bardir,$suffix) = fileparse($bar,qr/\.bar/);
  $bound = $bardir.$basename.".bound".$suffix;
  $unbound = $bardir.$basename.".unbound".$suffix;
}
else{
  $bound = $bar.".bound";
  $unbound = $bar.".unbound";
}
open $fhu, ">", $unbound or die "Cannot open filehandle for unbound $!";
open $fhb, ">", $bound or die "Cannot open filehandle for bound $!";

parse_barfile($bar);
consistify(\%lmins_u,$fhu);
consistify(\%lmins_b,$fhb);

close($fhu);
close($fhb);


sub consistify {
  my ($ref,$fh) = @_;
  my %basin = %$ref;
  my $start = 1; # needed to check for first lmin
  my $ii=1;
  foreach my $i (sort {$a <=> $b} keys %basin){
    my $b_father = $bar_o[$i-1][3];
    if($start == 1){
      printf($fh "     %s\n", $sequence);
      if ($i > 1){ # this subtree doesnt start with lmin 1
	$bar_o[$i-1][0] = $ii;
	$basin{$i} = $ii;
      }
      $start=0;
    }
    $bar_o[$i-1][0] = $ii; # assign new lmin number
    $basin{$i}=$ii;        # map old lmin number -> new lmin number
    if($b_father == 0){    # adjust father
      $bar_o[$i-1][3] = 0;
    }
    else{
      $bar_o[$i-1][3] = $basin{$b_father};
    }
    print_barline($bar_o[$i-1],$fh);
    $ii++;
  }

}

sub print_barline {
  my ($ref,$fh) = @_;
  my @b = @$ref;
  printf($fh "%4d %s %6.2f %4d %6.2f %s %12ld %8ld %10.6f %8ld %10.6f\n", $b[0],$b[1],$b[2],$b[3],$b[4],$b[5],$b[6],$b[7],$b[8],$b[9],$b[10]);
  return;
}

sub parse_barfile {
  my ($barfile) = @_;
  my $fh = new IO::File "< $barfile" or die "no file: $barfile found\n";
  my $firstline = <$fh>;
  $sequence = $1 if ($firstline =~ /(\S+)/);

  my $idx = 0;
  my ($mdx, $min_dist, $mstr);
  while(<$fh>) {
    my @line = split;
    push @bar_o, \@line;
    if ($line[1] =~ /\*/){
      $lmins_b{$line[0]}=1;
    }
    else {
      $lmins_u{$line[0]}=1;
    }
  }
  $fh->close;
  return;
}

__END__

=head1 NAME

crossrates.pl - A barriers-ligand post-processor

=head1 SYNOPSIS

crossrates.pl [-b|--bar I<FILE>] [options]

=head1 DESCRIPTION

This script is a post-processing utility for "barriers -M ligand"
calls.

It adjusts splits ligand bound (aka 'starred') and unbound (regular)
states into separate bar files and adjusts for ligand concentration in
rates matrices produced by barriers.

=head1 OPTIONS

=over

=item B<-b|--bar>

RNA energy landscape in Barriers format (.bar file)


=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut
