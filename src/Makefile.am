if BUILD_SECIS_EXT
SEDIR=SECIS
endif

SUBDIRS = $(SEDIR) PoHo

bin_PROGRAMS = barriers

barriers_SOURCES = \
    main.c \
    hash_util.c \
    barriers.c \
    compress.c \
    stapel.c \
    ringlist.c \
    moves.c \
    trees.c \
    treeplot.c \
    utils.c \
    simple_set.c \
    cmdline.c \
    rnamoves.c

noinst_HEADERS = \
    barrier_types.h \
    barriers.h \
    hash.h \
    hash_util.h \
    pair_mat.h \
    ringlist.h \
    stapel.h \
    tree_types.h \
    trees.h \
    utils.h \
    treeplot.h \
    simple_set.h \
    compress.h \
    cmdline.h \
    rnamoves.h


#  uncomment the following if barriers requires the math library
if BUILD_SECIS_EXT
barriers_CPPFLAGS = \
    -I$(top_srcdir)/src/SECIS

barriers_LDADD= \
    -L$(top_builddir)/src/SECIS -lSECIS \
    $(RNALIB_LIBS) \
    -lm \
    -lstdc++
else
barriers_LDADD= \
    -lm
endif

EXTRA_DIST = \
    barriers.ggo


main.c: cmdline.h

cmdline.h cmdline.c: barriers.ggo
	gengetopt -i barriers.ggo

# remove all generated man pages
maintainer-clean-local:
	rm -rf cmdline.h cmdline.c
