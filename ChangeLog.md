# Changelog

Revision history for barriers

### 1.8.1   2019-11-21    Gregor Entzian, Ronny Lorenz    <rna@tbi.univie.ac.at>
  * Fix connected component option
  * Update connected component description
  * Fix various mixed usages of signed/unsigned types
  * Replace hash table implementation with one that uses separate chaining with linked lists instead of simple linear probing
  * Refrain from installing `get_hashbits.pl` since this is only required at compile-time
  * Refactor installation process of Perl 5 scripts and Perl 5 module


### 1.8.0   2019-10-18    Gregor Entzian, Ronny Lorenz    <rna@tbi.univie.ac.at>
  * Fix command line argument parsing
  * Fix potentially running into infinite loops
  * Change data type for minima indices to unsigned long
  * Fix connected component output and help message
  * Add command line parameters to specify rate file names
  * Abort if number of input lines exceeds 2/3 of hash size
  * Re-integrate command line parameter to specify mapped structure output
  * Fix several memory leaks
  * Always print graph and move set options to STDERR unless --quiet option is provided
  * Use Markdown files for ChangeLog and README
  * Add alternative `RNA2` graph and move set implementation by Felix Kuehnl


### 1.7.0   2017-11-26    Michael T. Wolfinger    <mtw@tbi.univie.ac.at>
  * Disable Shift moves by default
  * Added noLP-rate functionality
  * Added ligand binding functionality
  * Added /scripts/crossrates.pl for post-processing ligand data
  * Allow up to 34 HASHBITS in conigure.ac (default 27)


### 0.01    2001-03-08    Ivo L. Hofacker, Christoph Flamm, Peter F. Stadler    <ivo@tbi.univie.ac.at>
  * barriers: initial version.
