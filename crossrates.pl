#!/usr/bin/perl
# -*-CPerl-*-
# Last changed Time-stamp: <2017-10-20 15:48:54 mtw>

use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use File::Basename;
use strict;
use warnings;

my ($sequence,$fnunbound,$fnbound,$fhu,$fhb,$fhra,$fhrb);
my ($ratesr,$dim,$fnrbr,$fhrbr);
my $bar=undef;
my $rates="rates.out";
my $binrates="rates.bin";
my %lmins_b = (); # indices of ligand-bound lmins
my %lmins_u = (); # indices of unbound lmins
my @bar_o   = (); # AoH holding the original bar file
my $T=37.;
my $K0=-273.15;
my $DuplexInit37=410;
my $DuplexInitdH=360;
my $DuplexInit = 0.;
my $Conc = 1.;
my $BindingBonus = 0.;


Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions(
					   "b|bar=s"    => \$bar,
					   "rates=s"    => \$rates,
					   "binrates=s" => \$binrates,
					   "T|temp=f"   => \$T,
					   "C|conc=f"   => \$Conc,
					   "E|energy=f" => \$BindingBonus,
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
  my ($basename,$dir,$suffix) = fileparse($bar,qr/\.bar/);
  $fnbound = $dir.$basename.".bound".$suffix;
  $fnunbound = $dir.$basename.".unbound".$suffix;
}
else{
  $fnbound = $bar.".bound";
  $fnunbound = $bar.".unbound";
}
if($binrates =~ m/\.bin$/){
  my ($basename,$dir,$suffix) = fileparse($binrates,qr/\.bin/);
  $fnrbr = $dir.$basename.".r".$suffix;
}
else{
  $fnrbr = $binrates.".r";
}

scale_duplexInitiationEnergy($T);
print "DuplexInit: $DuplexInit\n";
my $kT = 0.00198717*($K0+$T);
my $Beta = 1/$kT;

open $fhu, ">", $fnunbound
  or die "Cannot open filehandle for unbound bar file: $!\n";
open $fhb, ">", $fnbound
  or die "Cannot open filehandle for bound bar file: $!\n";
$dim = parse_barfile($bar);
print "dim is $dim\n";
consistify(\%lmins_u,$fhu);
consistify(\%lmins_b,$fhb);
close($fhu);
close($fhb);
#print Dumper(\%lmins_u);
#print Dumper(\%lmins_b);
open $fhra, "<", $rates
  or die "Cannot open filehandle for ASCII rates file: $!\n";
close($fhra);

open $fhrb, "<:raw", $binrates
  or die "Cannot open filehandle for reading binary rates file: $!\n";
open $fhrbr, ">:raw", $fnrbr
  or die "Cannot open filehandle for writing binary rates file: $!\n";
$rates = parse_binrates($fhrb);
$ratesr = reorder_matrix($rates,\%lmins_u,\%lmins_b,$dim);
#dump_matrix($ratesr,$dim);
write_binrates($ratesr,$fhrbr,$dim);
close($fhrb);
close($fhrbr);

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
  return \@m;
}

sub write_binrates {
  my ($r,$fh,$d) = @_;
  print $fh pack("i", $d);
  for (my $i=0;$i<$d;$i++){
    for (my $j=0;$j<$d;$j++){
      print $fh pack("d", $$r[$d*$i+$j]);
    }
  }
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
  
  print "LU:\n";
  print Dumper($lu);
  print "LB:\n";
  print Dumper($lb);

  # merge lookup lmin hashes
  foreach my $key (keys %$lu){ # unbound states
    if($$lu{$key}>$maxi){$maxi =  $$lu{$key}}
    $lm{$key} = $$lu{$key};
  }
  print "maxi=$maxi\n";
  foreach my $key (keys %$lb){ # bound states
    die "inconsistency in lmin hash: lmin $key already assigned\n"
      if (exists $lm{$key});
    $lm{$key} = $$lb{$key}+$maxi;
  }

  #print Dumper(\%lm);
  %lmr = reverse %lm;
  #print Dumper(\%lmr);
  for($i=1;$i<=$dim;$i++){
    $oi=$lmr{$i};
    $oi--;
    for($j=1;$j<=$dim;$j++){
      $oj=$lmr{$j}-1;
      $n[$dim*($i-1)+($j-1)]=$$r[$dim*$oi+$oj];
    }
#    print "\n";
  }

  my $a = adjust_crossterms(\@n,$maxi,$dim);
#  dump_matrix($a,$dim);

  
  return $a;
}

sub  adjust_crossterms{
  my ($mx,$max,$dim) = @_;
  # dump_matrix($mx,$dim);
  for (my $i=0;$i<$max;$i++){ # ON rate
    for(my $j=$max;$j<$dim;$j++){
      # $$mx[$dim*$i+$j] *= $Conc*$DuplexInit;
      $$mx[$dim*$j+$i] *= $Conc*$DuplexInit;
    }
  }

   for (my $i=$max;$i<$dim;$i++){ # OFF rate
    for(my $j=0;$j<$max;$j++){
      # $$mx[$dim*$i+$j] *= $DuplexInit * exp(-$Beta*$BindingBonus);
      $$mx[$dim*$j+$i] *= $DuplexInit * exp(-$Beta*$BindingBonus);
    }
  }
  # dump_matrix($mx,$dim);
  return $mx;
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


sub scale_duplexInitiationEnergy {
  my $T = shift;
  my $tempf  = (($T + $K0) / (37 + $K0));
  $DuplexInit = $DuplexInitdH - ($DuplexInitdH - $DuplexInit37) * $tempf;
  $DuplexInit /= 100.;
}

__END__

=head1 NAME

crossrates.pl - A barriers-ligand post-processor

=head1 SYNOPSIS

crossrates.pl [-b|--bar I<FILE>] [options]

=head1 DESCRIPTION

This script is a post-processing utility for "barriers -M ligand"
calls.

It deconvolutes ligand bound (aka 'starred') and unbound (regular)
states into separate bar files, reorganizes barrier rates matrices and
adjusts them for different ligand binding energies, concentrations.

=head1 OPTIONS

=over

=item B<-b|--bar>

RNA energy landscape in Barriers format (.bar file)

=item B<-T|--temp>

Simulation temperature in Celsius (default: 37.0)

=item B<--binrates>

Binary rates file (default: 'rates.bin')

=item B<--rates>

ASCII rates file (default: 'rates.out')

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut
