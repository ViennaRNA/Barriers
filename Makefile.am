AUTOMAKE_OPTIONS = foreign

if BUILD_SECIS_EXT
SEDIR=SECIS
endif

SUBDIRS = src $(SEDIR) PoHo man tests

bin_SCRIPTS = \
    scripts/crossrates.pl \
    scripts/get_hashbits.pl


EXTRA_DIST = \
    README.md \
    ChangeLog.md \
    barriers.lsm.in \
    VERSION.txt \
    scripts/crossrates.pl \
    scripts/get_hashbits.pl
