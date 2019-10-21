AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_builddir)/src/

# Unit Tests for RNA moveset
# Compile only if Check Unit Testing Framework found
if CHECK_FOUND
TESTS = check_rnamoves
check_PROGRAMS = check_rnamoves
check_rnamoves_SOURCES = check_rnamoves.c $(top_srcdir)/src/stapel.c $(top_srcdir)/src/utils.c $(top_srcdir)/src/ringlist.c # rnamoves.c 
check_rnamoves_CFLAGS = @CHECK_CFLAGS@
check_rnamoves_LDADD = @CHECK_LIBS@
endif
