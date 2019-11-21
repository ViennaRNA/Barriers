#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

# ~~~~~~~~~~~~~~~~ #
# GLOBAL VARIABLES #
# ~~~~~~~~~~~~~~~~ #
my $HB;  # HashBits
my $RAM; # RAM Size of Machine
my $HE;  # Bytes per Hash Entry
my $SL;  # Maximum Stringlength
my $CL;  # Compression Level
my $BS;  # Bytes per compressed String

$HE  = 48;
$CL  = 7;

&usage() unless GetOptions(
  "compression=i" => \$CL,
  "hashentry=i"   => \$HE,
  "ram=i"         => \$RAM,
  "maxlen=i"      => \$SL,
  "hashbits=i"    => \$HB,
  "v");

usage() unless $RAM && $SL;

$BS = roundup($SL/$CL)+1;

$HB = int(log2(($RAM*1024)/($HE+$BS))) unless $HB;

print STDERR "With $HB HASHBITS you have:\n".
"*) a hash that can store ".(2**$HB)." structures of length <= ".(($BS-1)*$CL)."\n",
"*) ".(($RAM*1024 - (2**$HB*($HE+$BS)))/1024)." kB of RAM left for the system.\n";

print "./configure --with-hash-bits=".$HB." CFLAGS=\"-mcmodel=large\"\n";

sub usage {
  die "\n$0 [-ram r] [-maxlen m]\n\n",
  "Required parameters:\n",
  "-ram r    \tRAM of the machine in kB\n",
  "-maxlen m \tmaximal sequence length\n\n",
  "Optional parameters:\n",
  "-compression\tString Compression Level of barriers (Default: $CL)\n",
  "-hashentry  \tAdditional bytes per Hashentry reserved by the System (Default: $HE)\n",
  "-hashbits   \tCompute RAM requirements for given number of Hashbits\n",
  "\n~~ HASHBITS = log2(ram(bytes)/($HE+((maxlen/$CL)+1))) ~~\n\n";
}

sub log2 {
  my $n = shift;
  return log($n)/log(2);
}

sub roundup {
  my $n = shift;
  return (($n == int($n)) ? $n : int($n+1));
}
