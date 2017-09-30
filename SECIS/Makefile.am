# $Id: Makefile.am,v 1.1 2006/07/18 14:08:49 xtof Exp $
# Last changed Time-stamp: <2017-09-30 18:16:49 xtof>
AM_CPPFLAGS = -I@ac_VRNA_includes@
LDADD = -L@ac_VRNA_lib@ -lRNA -lstdc++
lib_LIBRARIES = libSECIS.a
libSECIS_a_SOURCES = secis_basics.cpp \
                     secis_constraints.cpp \
	             secis_generate.cpp \
                     secis_neighbors.cpp \
	             secis_struct.cpp
pkginclude_HEADERS = secis_basics.h \
	             secis_constraints.h \
	             secis_generate.h \
                     secis_neighbors.h \
	             secis_struct.h
bin_PROGRAMS = secisFindSeqs
secisFindSeqs_SOURCES = secis_find_seq.cpp
secisFindSeqs_LDADD = -L. -lSECIS $(LDADD)
EXTRA_DIST = Examples/m9s30.bar \
             Examples/m9s30_bar.in \
             Examples/m9s30_tree.ps \
             Examples/README \
             Examples/secis.txt \
             Examples/sequences_m9_s30 \
             Examples/to_bar_in.pl

