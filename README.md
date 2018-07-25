`barriers` - Compute local minima and energy barriers of a landscape

## INSTALLATION

Typically, `barriers` in installed like this:

```[bash]
 ./configure
 ./make
 ./make install
```

The `barriers` program uses a large static array for the hash table at the
heart of the method. The size of this hash and therefore the maximum size
of the landscape you can analyse is set at compile time using the
`--with-hash-bits` option of the configure script, e.g.:
```[bash]
 ./configure --with-hash-bits=27
```
will create a hash of 2^27 entries.

If you you need a hash size >2^27 you may have to set the `-mcmodel=medium`
or `-mcmodel=large` option in `gcc`, e.g.:
```[bash]
 ./configure --with-hash-bits=29 CFLAGS='-mcmodel=medium'
```

You can use the script `scripts/get_hashbits.pl` to optimize the hashsize
given your RAM size and sequence length.

## PUBLICATIONS

For an explanation of the algorithm, please see 

Ch. Flamm, I.L. Hofacker, P.F. Stadler, M.T. Wolfinger.
Barrier trees of degenerate landscapes. 
Z. Phys. Chem., 216:155–173, 2002. (doi)[https://doi.org/10.1524/zpch.2002.216.2.155]

Ligand support has been added with version 1.7.0, see

M.T. Wolfinger, Ch. Flamm, I.L. Hofacker
Efficient computation of cotranscriptional RNA-ligand interaction dynamics
Methods, 2018 (doi)[https://doi.org/10.1016/j.ymeth.2018.04.036]


## COPYRIGHT AND LICENCE

Copyright (C) 2001 Ivo Hofacker, Christoph Flamm, Peter Stadler, Michael
T. Wolfinger

`barriers` is free software; you can redistribute it and/or modify it under
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
