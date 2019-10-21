[![GitHub release](https://img.shields.io/github/release/ViennaRNA/Barriers.svg)](https://www.tbi.univie.ac.at/RNA/Barriers/#download)
[![Build Status](https://travis-ci.org/ViennaRNA/Barriers.svg?branch=master)](https://travis-ci.org/ViennaRNA/Barriers)
[![Github All Releases](https://img.shields.io/github/downloads/ViennaRNA/Barriers/total.svg)](https://github.com/ViennaRNA/Barriers/releases)
[![Conda](https://img.shields.io/conda/v/bioconda/barriers.svg)](https://anaconda.org/bioconda/barriers)
[![AUR](https://img.shields.io/aur/version/barriers.svg)](https://aur.archlinux.org/packages/barriers/)

# Barriers

Compute local minima and energy barriers of a landscape

----

## INSTALLATION

Typically, barriers in installed like this:

```
./configure
./make
./make install
```

The barriers program uses a large static array for the hash table at the
heart of the method. The size of this hash and therefore the maximum size
of the landscape you can analyse is set at compile time using the
`--with-hash-bits` option of the configure script, e.g.:

```
./configure --with-hash-bits=27
```

will create a hash of 2^27 entries.

If you you need a hash size >2^27 you may have to set the `-mcmodel=medium`
or `-mcmodel=large` option in gcc, e.g.:

```
./configure --with-hash-bits=29 CFLAGS='-mcmodel=medium'
```

You can use the script 'scripts/get_hashbits.pl' to optimize the hashsize
given your RAM size and sequence length.

----

## IMPORTANT NOTES

With version 1.7.0, the default move set (`-M`) for RNA graphs changed from
Shift moves to just insertion/removal of single base pairs! To reproduce
results obtained with barriers prior to 1.7.0, you need to explicitly
activate the Shift move option, e.g.:

```
barriers -M Shift
```

----

## PUBLICATIONS

For an explanation of the algorithm, please see 

Ch. Flamm, I.L. Hofacker, P.F. Stadler, M.T. Wolfinger.
Barrier trees of degenerate landscapes. 
Z. Phys. Chem., 216:155--173, 2002. (doi:10.1524/zpch.2002.216.2.155)

Ligand support has been added with version 1.7.0, see

M.T. Wolfinger, Ch. Flamm, I.L. Hofacker
Efficient computation of cotranscriptional RNA-ligand interaction dynamics
Methods, 2018 (doi: 10.1016/j.ymeth.2018.04.036)

----

## COPYRIGHT AND LICENCE

Copyright (C) 2001 Ivo Hofacker, Christoph Flamm, Peter Stadler, Michael
T. Wolfinger

barriers is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc., 59
Temple Place, Suite 330, Boston, MA 02111-1307 USA

Comments are welcome <ivo@tbi.univie.ac.at>
