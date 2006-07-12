if BUILD_SECIS_EXT
SEDIR=SECIS
endif

SUBDIRS = $(SEDIR) .

bin_PROGRAMS=barriers
barriers_SOURCES=main.c hash_util.c barriers.c compress.c stapel.c \
        ringlist.c moves.c trees.c treeplot.c utils.c simple_set.c \
	getopt.c getopt1.c

noinst_HEADERS = barrier_types.h barriers.h hash.h hash_util.h pair_mat.h \
        ringlist.h stapel.h tree_types.h trees.h utils.h treeplot.h \
	getopt.h simple_set.h compress.h

#  uncomment the following if barriers requires the math library
if BUILD_SECIS_EXT
barriers_LDADD= -LSECIS -lSECIS -L@ac_VRNA_lib@ -lRNA -lm
else
barriers_LDADD= -lm 
endif

EXTRA_DIST=barriers.lsm.in barriers.spec.in barriers.texinfo barriers.1

#  if you write a self-test script named `chk', uncomment the
#  following and add `chk' to the EXTRA_DIST list
#TESTS=chk

#  build and install the .info pages
# info_TEXINFOS = barriers.texinfo
# barriers_TEXINFOS = gpl.texinfo

#  install the man pages
man_MANS=barriers.1
