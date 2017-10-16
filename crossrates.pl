#!/usr/bin/perl
# -*-CPerl-*-
# Last changed Time-stamp: <2017-10-16 17:44:00 mtw>

use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use File::Basename;
use strict;
use warnings;

my ($sequence,$unbound,$bound,$fhu,$fhb,$fhra,$fhrb);
my ($basename,$bardir,$suffix,$ratesr,$dim);
my $bar=undef;
my $rates="rates.out";
my $binrates="rates.bin";
my %lmins_b = (); # indices of ligand-bound lmins
my %lmins_u = (); # indices of unbound lmins
my @bar_o   = (); # AoH holding the original bar file

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions(
					   "b|bar=s"    => \$bar,
					   "rates=s"    => \$rates,
					   "binrates=s" => \$binrates,
					   "man"        => sub{pod2usage(-verbose => 2)},
                                           "h|help"     => sub{pod2usage(1)}
					  );
unless (-f $rates) {
  warn "Could not find barriers ASCII rates file '$rates'";
  pod2usage(-verbose => 0);
}
unless (-f $binrates) {
  warn "Could not find barriers binary rates file '$binrates'";
  pod2usage(-verbose => 0);
}
unless (-f $bar){
  warn "Could not find barriers .bar file '$bar'";
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
open $fhu, ">", $unbound
  or die "Cannot open filehandle for unbound: $!\n";
open $fhb, ">", $bound
  or die "Cannot open filehandle for bound: $!\n";
$dim = parse_barfile($bar);
consistify(\%lmins_u,$fhu);
consistify(\%lmins_b,$fhb);
close($fhu);
close($fhb);
#print Dumper(\%lmins_u);
#print Dumper(\%lmins_b);
open $fhra, "<", $rates
  or die "Cannot open filehandle for ASCII rates file: $!\n";
open $fhrb, "<:raw", $binrates
  or die "Cannot open filehandle for binary rates file: $!\n";
$rates = parse_binrates($fhrb);
$ratesr = reorder_matrix($rates,\%lmins_u,\%lmins_b,$dim);
close($fhra);
close($fhrb);

sub consistify {
  my ($B,$fh) = @_;
  my $start = 1;
  my $ii = 1;

  foreach my $i (sort {$a <=> $b} keys %$B){
    my $b_father = $bar_o[$i-1][3];
    if($start == 1){ # for first lmin only
      printf($fh "     %s\n", $sequence);
      if ($i > 1){   # this subtree doesnt start with lmin 1
	$bar_o[$i-1][0] = $ii;
	$$B{$i} = $ii;
      }
      $start=0;
    }
    $bar_o[$i-1][0] = $ii; # assign new lmin number
    $$B{$i}=$ii;           # map old lmin number -> new lmin number
    if($b_father == 0){    # adjust father
      $bar_o[$i-1][3] = 0;
    }
    else{
      $bar_o[$i-1][3] = $$B{$b_father};
    }
    print_barline($bar_o[$i-1],$fh);
    $ii++;
  }
}

sub print_barline {
  my ($ref,$fh) = @_;
  my @b = @$ref;
  printf($fh "%4d %s %6.2f %4d %6.2f %s %12ld %8ld %10.6f %8ld %10.6f\n",
	 $b[0],$b[1],$b[2],$b[3],$b[4],$b[5],$b[6],$b[7],$b[8],$b[9],$b[10]);
  return;
}

sub parse_barfile {
  my ($barfile) = @_;
  my $fh = new IO::File "< $barfile" or die "no file: $barfile found\n";
  my $firstline = <$fh>;
  my $i = 0;
  $sequence = $1 if ($firstline =~ /(\S+)/);

  my $idx = 0;
  my ($mdx, $min_dist, $mstr);
  while(<$fh>) {
    $i++;
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
  return $i;
}

sub parse_binrates {
  my $fh = shift;
  my @m=();
  my $d;

  undef $/;
  @m = unpack("id*", <$fh>);
  $d = shift @m;
  print "$d\n";
  for my $row (0..$d-1){
    for my $col (0..$d-1){
   #   printf (STDOUT "%10.4g ", $m[$d*$row+$col]);
    }
#    printf (STDOUT "\n");
  }
  return \@m;
}

sub reorder_matrix {
  my ($r,$lu,$lb,$d) = @_;
  my $maxi=1;
  my ($i,$j,$oi,$oj);
  my @n = ();
  my %lm = (); # merged lmin map old id -> new id
  my %lmr = ();

  # initialize new, reordered rates matrix
  for(my $i=0;$i<$dim;$i++){
    for(my $j=0;$j<$dim;$j++){
      $n[$dim*$i+$j]=0.;
    }
  }

  # merge lookup lmin hashes
  foreach my $key (keys %$lu){ # unbound states
    if($$lu{$key}>$maxi){$maxi =  $$lu{$key}}
    $lm{$key} = $$lu{$key};
  }
  foreach my $key (keys %$lb){ # bound states
    die "inconsistency if lmin hash: lmin $key already assigned\n"
      if (exists $lm{$key});
    $lm{$key} = $$lb{$key}+$maxi;
  }
  #print Dumper(\%lm);
  %lmr = reverse %lm;
  #print Dumper(\%lmr);
  for($i=1;$i<=$dim;$i++){
    $oi=$lmr{$i};
    #print "oi=$oi\n";
    $oi--;
    for($j=1;$j<=$dim;$j++){
      $oj=$lmr{$j}-1;
      $n[$dim*($i-1)+($j-1)]=$$r[$dim*$oi+$oj];
    }
#    print "\n";
  }
  dump_matrix(\@n,$dim);
}

sub dump_matrix {
  my ($mx,$dim) = @_;
  for(my $i=0;$i<$dim;$i++){
    for(my $j=0;$j<$dim;$j++){
      printf STDOUT "%10.4g ",$$mx[$dim*$i+$j];
    }
    print STDOUT "\n";
  }
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

=item B<--binrates>

Binary rates file (default: 'rates.bin')

=item B<--rates>

ASCII rates file (default: 'rates.out')

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut
