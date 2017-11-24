#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2017-11-23 23:06:48 mtw>

use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use File::Basename;
use YAML::XS qw(LoadFile DumpFile);
use strict;
use warnings;

my ($fnbu,$fnbb,$fnbr,$fhbu,$fhbb,$fhbr,$fhra,$fhrb,$fnrbr,$fhrbr,$fnm);
my ($sequence,$ratesr,$ratesra,$dim,$dimr,$maxb,$rates);
my $bar=undef;
my $asciirates="rates.out";
my $binrates="rates.bin";
my %lminB    = (); # indices of ligand-bound lmins
my %lminU    = (); # indices of unbound lmins
my %lminMap  = (); # map old => new lmin numbers
my %lminMapR = (); # map new => old lmin numbers
my @bar_o    = (); # AoH holding the original bar file
my @bar_b    = (); # AoH holding the original bar file (backup copy)
my $T=37.;
my $K0=-273.15;
my $DuplexInit37=410;
my $DuplexInitdH=360;
my $DuplexInit=0.;
my $Conc=1.;
my $BindingBonus=0.;
my $ymlminmap=undef;
my $has_map=0;
my $Aa=600;  # association/dissociation rate

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1)
  unless GetOptions(
		    "a|assrate=f" => \$Aa,
		    "b|bar=s"     => \$bar,
		    # "rates=s"     => \$asciirates,
		    "binrates=s"  => \$binrates,
		    "T|temp=f"    => \$T,
		    "C|conc=f"    => \$Conc,
		    "E|energy=f"  => \$BindingBonus,
		    "map=s"       => \$ymlminmap,
		    "man"         => sub{pod2usage(-verbose => 2)},
		    "h|help"      => sub{pod2usage(1)}
		   );

if (defined $ymlminmap){    # don't process bar file
  unless (-f $ymlminmap) {  # read lmin map from YAML instead
    warn "Could not find YAML lmin map file '$ymlminmap'";
    pod2usage(-verbose => 0);
  }
  $has_map=1;
  my $lminMapRef = LoadFile($ymlminmap);
#  print Dumper($lminMapRef);
  my $xref = $$lminMapRef{map};
  %lminMap = %$xref;
  %lminMapR = reverse %lminMap;
  $maxb = %$lminMapRef{unbound};
  $dim = %$lminMapRef{dim};
#  print Dumper(\%lminMap);
#  print Dumper($unbound);
#  die;
}
#unless (-f $asciirates) {
#  warn "Could not find barriers ASCII rates file '$rates'";
#  warn " ... continuing with binary rates";
  # pod2usage(-verbose => 0);
#}
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
  $fnbb = $dir.$basename.".b".$suffix;
  $fnbu = $dir.$basename.".u".$suffix;
  $fnbr = $dir.$basename.".r".$suffix;
  $fnm  = $dir.$basename.".r.yaml";
}
else{
  $fnbb = $bar.".b";
  $fnbu = $bar.".u";
  $fnbr = $bar.".r";
  $fnm  = $bar.".r.yaml";
}
if($binrates =~ m/\.bin$/){
  my ($basename,$dir,$suffix) = fileparse($binrates,qr/\.bin/);
  $fnrbr = $dir.$basename.".r".$suffix;
}
else{$fnrbr = $binrates.".r"}

scale_duplexInitiationEnergy($T);
#print "* crossrates.pl: DuplexInit=$DuplexInit\n";
my $kT = 0.00198717*($K0+$T);
my $Beta = 1/$kT;

unless ($has_map == 1){ # if YAML lmin map was not provided
  open $fhbu, ">", $fnbu
    or die "Cannot open filehandle for unbound bar file: $!\n";
  open $fhbb, ">", $fnbb
    or die "Cannot open filehandle for bound bar file: $!\n";
  open $fhbr, ">", $fnbr
    or die "Cannot open filehandle for reordered bar file: $!\n";
  $dim = parse_barfile($bar);
  consistify(\%lminU,$fhbu);
  consistify(\%lminB,$fhbb);
  $maxb = lmins_old2new(\%lminU,\%lminB,$dim);
  write_reordered_barfile(\%lminMap,\%lminMapR,$fhbr);
  close($fhbu);
  close($fhbb);
  close($fhbr);
}

#open $fhra, "<", $rates
#  or die "Cannot open filehandle for ASCII rates file: $!\n";
#close($fhra);
open $fhrb, "<:raw", $binrates
  or die "Cannot open filehandle for reading binary rates file: $!\n";
open $fhrbr, ">:raw", $fnrbr
  or die "Cannot open filehandle for writing binary rates file: $!\n";
($dimr,$rates) = parse_binrates($fhrb);
die "dimension mismatch: $dim from YAML file != $dimr from rates file $!\n"
  unless ($dim == $dimr);
$ratesr = reorder_matrix($rates,$dim); # re-order matrix, unbound states first
$ratesra = adjust_crossterms($ratesr,$maxb,$dim);
#dump_matrix($ratesra,$dim);
write_binrates($ratesra,$fhrbr,$dim);
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

sub write_reordered_barfile {
  my ($B,$BR,$fh) = @_;
  printf($fh "     %s\n", $sequence);
  foreach my $i (sort {$a <=> $b} keys %$BR){
    my $olmin = $$BR{$i};
    $bar_b[$olmin-1][0] = $i;
    my $ofather =  $bar_b[$olmin-1][3];
    if (defined $$B{$ofather}){
      $bar_b[$olmin-1][3] = $$B{$ofather};
    }
    else{
      $bar_b[$olmin-1][3] = 0;
    }
    print_barline($bar_b[$olmin-1],$fh);
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
      $lminB{$line[0]}=1;
    }
    else {
      $lminU{$line[0]}=1;
    }
  }
  $fh->close;
  @bar_b = map { [@$_] } @bar_o; # deep copy bar_o
  return $i;
}

sub parse_binrates {
  my $fh = shift;
  my @m=();
  my $d;

  undef $/;
  @m = unpack("id*", <$fh>);
  $d = shift @m;
 # print "$d\n";
  return ($d, \@m);
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

sub lmins_old2new {
  my ($lu,$lb,$d) = @_;
  my $maxi=1;

  # merge lookup lmin hashes
  foreach my $key (keys %$lu){ # unbound states
    if($$lu{$key}>$maxi){$maxi = $$lu{$key}}
    $lminMap{$key} = $$lu{$key};
  }
  foreach my $key (keys %$lb){ # bound states
    die "inconsistency in lmin hash: lmin $key already assigned\n"
      if (exists $lminMap{$key});
    $lminMap{$key} = $$lb{$key}+$maxi;
  }

  %lminMapR = reverse %lminMap;

  my %lminH = (unbound => $maxi,
	       dim     => $d,
	       map     => \%lminMap);

  DumpFile($fnm, \%lminH);

  # print "LU:\n"; print Dumper($lu); print "LB:\n"; print Dumper($lb);
  # print "LM:\n"; print Dumper(\%lminMap);
  # print "LMR:\n";  print Dumper(\%lminMapR);
  return $maxi;
}

sub reorder_matrix {
  my ($r,$d) = @_;
  my ($i,$j,$oi,$oj);
  my @n = ();

  for($i=0;$i<$dim;$i++){
    for($j=0;$j<$dim;$j++){
      $n[$dim*$i+$j]=0.;  # initialize new, reordered rates matrix
    }
  }

  for($i=1;$i<=$dim;$i++){
    $oi=$lminMapR{$i};
    for($j=1;$j<=$dim;$j++){
      $oj=$lminMapR{$j};
      $n[$dim*($i-1)+($j-1)]=$$r[$dim*($oi-1)+($oj-1)];
    }
  }
  return \@n;
}

sub  adjust_crossterms{
  my ($mx,$max,$dim) = @_;

  # dump_matrix($mx,$dim);
  for (my $i=0;$i<$max;$i++){ # ON rate
    for(my $j=$max;$j<$dim;$j++){
      $$mx[$dim*$j+$i] *= $Conc*$Aa;
    }
  }

   for (my $i=$max;$i<$dim;$i++){ # OFF rate
    for(my $j=0;$j<$max;$j++){
      $$mx[$dim*$j+$i] *= $Aa * exp(-$Beta*$BindingBonus);
    }
  }
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

crossrates.pl [-b|--bar I<FILE>] [--binrates I<FILE>] [options]

=head1 DESCRIPTION

This script is a post-processing utility for "barriers -M ligand"
calls. It expects a bar file provided via the B<'--bar'> command line
option as well as a I<binary> rates matrix provided via the
B<'--binrates'> option.

It deconvolutes ligand bound (aka 'starred') and unbound (regular)
states into separate bar files, and produces a reorganized bar file
where unbound states are listed first, followed by all bound
states. It further reorganizes barriers rates matrices accordingly and
adjusts the cross-terms (i.e. those holding k_on and k_off rates)
for different ligand binding energies and concentrations. This
reorganization is done primarily to facilitate post-processing by
e.g. treekin or BarMap by grouping (and thereby separating) bound and
unbound states.

This script also creates a YAML file that contains mapping information
on lmin reorganization. When called with the B<--map> option, such a
YAML file is expected and no bar / rates file reorganization is
performed. This is espeically useful for adjusting the rates matrix'
cross terms for different ligand concentration and binding bonus
energies I<without> repleatedly performing the lmin reorganization.

=head1 OPTIONS

=over

=item B<-a|--assrate>

Association/dissociation rate (default: 600)

=item B<-b|--bar>

RNA energy landscape in Barriers format (.bar file)

=item B<-C|--conc>

Ligand concentration (default: 1.0)

=item B<-E|--energy>

Ligand binding (stabilization) energy in kcal/mol (default: 0.0)

=item B<-T|--temp>

Simulation temperature in Celsius (default: 37.0)

=item B<--binrates>

Binary rates file (default: 'rates.bin')

=item B<--map>

Provide a (previously generated) YAML file containing the lmin
reorganization information. If this is given, no bar / rates file
reorganization is performed.


=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut
