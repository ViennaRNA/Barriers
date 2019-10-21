AUTOMAKE_OPTIONS = foreign

AM_DISTCHECK_CONFIGURE_FLAGS = \
    --with-secis \
    --disable-silent-rules

SUBDIRS = src man tests

bin_SCRIPTS = \
    perl/crossrates.pl \
    perl/get_hashbits.pl


EXTRA_DIST = \
    README.md \
    ChangeLog.md \
    VERSION.txt \
    perl/crossrates.pl \
    perl/get_hashbits.pl
