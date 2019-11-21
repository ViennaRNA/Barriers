
common_scripts = \
    crossrates.pl \
    saddle2dot.pl \
    saddle2gml.pl \
    treeplot.pl

if INSTALL_PERL

if INSTALL_PERL_MOD
perlmoddir = $(prefix)$(PERL_LIB_RELATIVE_INSTALL_DIR)

nobase_dist_perlmod_DATA = \
  RNA/barrier.pm

dist_bin_SCRIPTS = \
    $(common_scripts) \
    fix_bar.pl \
    saddle.pl
else
dist_bin_SCRIPTS = \
    $(common_scripts)

endif

EXTRA_DIST = \
    get_hashbits.pl

else # without any perl scripts or modules but we still distribute files

EXTRA_DIST = \
    $(common_scripts) \
    get_hashbits.pl

endif # INSTALL_PERL

