AUTOMAKE_OPTIONS = foreign

if BUILD_SECIS_EXT
SEDIR=SECIS
endif

SUBDIRS = src $(SEDIR) PoHo man tests

bin_SCRIPTS = \
    perl/crossrates.pl \
    perl/get_hashbits.pl


EXTRA_DIST = \
    README.md \
    ChangeLog.md \
    VERSION.txt \
    perl/crossrates.pl \
    perl/get_hashbits.pl
