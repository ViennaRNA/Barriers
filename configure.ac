dnl Process this file with autoconf to produce a configure script.

dnl define package name and version
AC_INIT([barriers],
	m4_esyscmd_s([cat VERSION.txt]),
	[{ivo,xtof,mtw}@tbi.univie.ac.at],[Barriers], [https://www.tbi.univie.ac.at/RNA/Barriers])

AM_INIT_AUTOMAKE

dnl create a config.h file (Automake will add -DHAVE_CONFIG_H)
AM_CONFIG_HEADER(config.h)

dnl get current date and time
ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)

AC_CANONICAL_HOST

dnl configure options
AM_WITH_DMALLOC
hashbits=24
AC_ARG_WITH(hash-bits,
[  --with-hash-bits=ARG    use hashtable with 2^ARG entries (<32, default 24)],
[
  hashbits=$withval
  if test $hashbits -ge 32; then
    AC_MSG_WARN([Hash table size too large, setting to default size=24])
    hashbits=24
  else
    if test $hashbits -ge 27; then
      CFLAGS="$CFLAGS -mcmodel=large"
    fi
  fi
],
)
AC_DEFINE_UNQUOTED(HASHBITS, $hashbits, [use hash with 2^HASHBITS entries])

AC_ARG_WITH(secis,
[  --with-secis          build SECIS element design extension],
[ with_secis="yes" ],
[ with_secis="no"]
)

dnl default install location for Vienna RNA Package
ac_VRNA_includes=/usr/local/ViennaRNA
ac_VRNA_lib=/usr/local/lib

AC_ARG_WITH(ViennaRNA-include,
  [ --with-ViennaRNA-include    where to search for ViennaRNA header files],
  ac_VRNA_includes="$withval",
  )

AC_ARG_WITH(ViennaRNA-lib,
  [ --with-ViennaRNA-lib        where to search for libRNA.a],
  ac_VRNA_lib="$withval",
  )

AC_SUBST(ac_VRNA_includes)
AC_SUBST(ac_VRNA_lib)

dnl Checks for programs.
AC_PROG_INSTALL
AC_PROG_CC
AC_PROG_CXX
AC_PROG_CPP
AC_PROG_RANLIB

AC_C_INLINE

dnl Checks for libraries.

dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS(string.h unistd.h)

dnl Checks for typedefs, structures, and compiler characteristics.

dnl Checks for library functions.
AC_CHECK_FUNCS(erand48)

dnl Conditionally build Makefile in SECIS subdirectory
have_secis_ext=0
if test "$with_secis" != "no"; then
   if test -f $srcdir/SECIS/Makefile.am; then
      have_secis_ext=1 
dnl      AC_OUTPUT(SECIS/Makefile)
   else
      AC_MSG_WARN([SECIS subdirectory does not exist])
      with_secis="no"
   fi
fi

AC_DEFINE_UNQUOTED(HAVE_SECIS_EXTENSION,
	 $have_secis_ext,
	 [build with SECIS design extension])	

AM_CONDITIONAL(BUILD_SECIS_EXT, test "$with_secis" != "no")

AC_OUTPUT([Makefile barriers.lsm])
AC_OUTPUT([PoHo/Makefile])
AC_OUTPUT([SECIS/Makefile])
AC_OUTPUT([packaging/barriers.spec packaging/PKGBUILD])


if test "$with_secis" != "no"; then
AC_MSG_NOTICE(
[
Configure successful with the following options:

  SECIS Extension:               ${with_secis:-yes}
  Dependency to ViennaRNA        ${with_secis:-yes}
  Expect ViennaRNA headers in:   $ac_VRNA_includes
  Expect ViennaRNA library in:   $ac_VRNA_lib

If the ViennaRNA headers and library are in a different location on
your system, please set the appropriate absolute paths with the
configure options --with-ViennaRNA-include and --with-ViennaRNA-lib
])
fi
