if BUILD_SECIS_EXT
SEDIR=SECIS
endif

SUBDIRS = $(SEDIR) PoHo .

bin_PROGRAMS=barriers
barriers_SOURCES=main.c hash_util.c barriers.c compress.c stapel.c \
        ringlist.c moves.c trees.c treeplot.c utils.c simple_set.c \
	cmdline.c rnamoves.c

noinst_HEADERS = barrier_types.h barriers.h hash.h hash_util.h pair_mat.h \
        ringlist.h stapel.h tree_types.h trees.h utils.h treeplot.h \
	simple_set.h compress.h cmdline.h rnamoves.h

# Unit Tests for RNA moveset
# Compile only if Check Unit Testing Framework found
if CHECK_FOUND
TESTS = check_rnamoves
check_PROGRAMS = check_rnamoves
check_rnamoves_SOURCES = check_rnamoves.c stapel.c utils.c ringlist.c
check_rnamoves_CFLAGS = @CHECK_CFLAGS@
check_rnamoves_LDADD = @CHECK_LIBS@
endif

bin_SCRIPTS = scripts/crossrates.pl scripts/get_hashbits.pl
# scripts/crossrates.pl: scripts/crossrates.pl
# scripts/get_hashbits.pl: scripts/get_hashbits.pl

#  uncomment the following if barriers requires the math library
if BUILD_SECIS_EXT
barriers_LDADD= -LSECIS -lSECIS -L@ac_VRNA_lib@ -lRNA -lm -lstdc++
else
barriers_LDADD= -lm 
endif

EXTRA_DIST=barriers.lsm.in barriers.texinfo barriers.1 barriers.ggo \
	VERSION.txt scripts/crossrates.pl scripts/get_hashbits.pl

#  if you write a self-test script named `chk', uncomment the
#  following and add `chk' to the EXTRA_DIST list
#TESTS=chk

#  build and install the .info pages
# info_TEXINFOS = barriers.texinfo
# barriers_TEXINFOS = gpl.texinfo

#  install the man pages
man_MANS=barriers.1

main.c: cmdline.h


cmdline.h cmdline.c: barriers.ggo
	gengetopt -i barriers.ggo
