bin_PROGRAMS=barriers
barriers_SOURCES=main.c hash_util.c barriers.c stapel.c \
        ringlist.c moves.c trees.c treeplot.c utils.c \
	barriers.1 system.h aclocal.m4 getopt.c getopt1.c getopt.h \
	barrier_types.h barriers.h hash.h hash_util.h pair_mat.h \
	ringlist.h stapel.h tree_types.h trees.h utils.h treeplot.h

#  uncomment the following if barriers requires the math library
barriers_LDADD=-lm

EXTRA_DIST=barriers.lsm.in barriers.spec.in barriers.texinfo

#  if you write a self-test script named `chk', uncomment the
#  following and add `chk' to the EXTRA_DIST list
#TESTS=chk

#  build and install the .info pages
# info_TEXINFOS = barriers.texinfo
# barriers_TEXINFOS = gpl.texinfo

#  install the man pages
man_MANS=barriers.1
