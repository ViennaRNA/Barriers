dnl Process this file with autoconf to produce a configure script.

dnl define package name and version
AC_INIT([barriers],
	m4_esyscmd_s([cat VERSION.txt]),
	[{ivo,xtof,mtw}@tbi.univie.ac.at],[Barriers], [https://www.tbi.univie.ac.at/RNA/Barriers])

AC_PREREQ([2.60])
AC_CONFIG_AUX_DIR(config)

AC_SUBST(VERSION)

AC_CONFIG_MACRO_DIR([m4])

AM_SILENT_RULES([yes])
AM_INIT_AUTOMAKE([1.11 -Wall -Werror])

dnl enable unit tests if 'Check' framework is installed
AC_CHECK_HEADERS(
    [check.h],
    [check_header_found=true],
    [
        check_header_found=false
        AC_MSG_WARN( [Header of 'Check' framework not found, skipping unit tests])
    ]
)
PKG_CHECK_MODULES(
    [CHECK],
    [check >= 0.9.4],
    [ check_found=true ],
    [
        AC_MSG_WARN( ['Check' framework not found, skipping unit tests])
        check_found=false
    ]
)
AM_CONDITIONAL(
    [CHECK_FOUND],
    [ test x$check_header_found = xtrue &&
      test x$check_found = xtrue ])

dnl create a config.h file (Automake will add -DHAVE_CONFIG_H)
AM_CONFIG_HEADER(config.h)

dnl get current date and time
ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)

AC_CANONICAL_HOST

dnl configure options
AM_WITH_DMALLOC
hashbits=27
AC_ARG_WITH(hash-bits,
[  --with-hash-bits=ARG    use hashtable with 2^ARG entries (<32, default 24)],
[
  hashbits=$withval
  if test $hashbits -ge 35; then
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

dnl Checks for programs.
AC_PROG_INSTALL
AC_PROG_CC
AC_PROG_CXX
AC_PROG_CPP
AC_PROG_RANLIB
AM_PROG_AR

AC_C_INLINE

dnl Checks for libraries.

dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS(string.h unistd.h)

dnl Checks for typedefs, structures, and compiler characteristics.

dnl Checks for library functions.
AC_CHECK_FUNCS(erand48)


dnl ---------------------------------------------------------------------------
dnl Checks for gengetopt
dnl ---------------------------------------------------------------------------

AC_ARG_VAR([GENGETOPT], [the 'gengetopt' program to generate command line argument parsers for executable programs])
AC_PATH_PROG([GENGETOPT], [gengetopt], [no])
AC_SUBST([GENGETOPT])
AM_CONDITIONAL(BARRIERS_AM_SWITCH_HAS_GENGETOPT, test "x$GENGETOPT" != "xno")


dnl ---------------------------------------------------------------------------
dnl Miscellaneous
dnl ---------------------------------------------------------------------------

AC_ARG_VAR([HELP2MAN], [the 'help2man' script to generate man pages from command line options of our executable programs])
AC_PATH_PROG([HELP2MAN], [help2man], [no])
AC_SUBST([HELP2MAN])
AM_CONDITIONAL(BARRIERS_AM_SWITCH_HAS_HELP2MAN, test "x$HELP2MAN" != "xno")


dnl Conditionally build Makefile in SECIS subdirectory
have_secis_ext=0
if test "$with_secis" != "no"; then
  #Search for the RNAlib2 and set variables for Makefiles if found.  _CFLAGS and _LIBS can be overriden with the modue Name from PKG_CHECK_MODULES
  PKG_CHECK_MODULES([RNALIB], [RNAlib2 >= 2.2.0],
                    [ AC_SUBST([RNALIB_CFLAGS])
                      AC_SUBST([RNALIB_LIBS])
                    ],
                    [])

   if test -f $srcdir/src/SECIS/Makefile.am; then
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

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile tests/Makefile])
AC_CONFIG_FILES([src/PoHo/Makefile])
AC_CONFIG_FILES([src/SECIS/Makefile])
AC_CONFIG_FILES([packaging/barriers.spec packaging/PKGBUILD packaging/barriers.lsm])
AC_CONFIG_FILES([man/Makefile])
AC_CONFIG_FILES([man/cmdlopt.sh],[chmod +x man/cmdlopt.sh])

AC_OUTPUT

# get directory paths

eval _bindir=$(eval printf "%s" $bindir)
eval _datadir=$(eval printf "%s" $datadir)
eval _mandir=$(eval printf "%s" $mandir)

m4_map_args([ AC_COLOR_RESULT_PACKAGE],
            [secis])
m4_map_args([ AC_COLOR_RESULT_SIMPLE],
            [_bindir],
            [_mandir],
            [_datadir])

# Notify the user

AC_STRING_APPEND_FORMAT_BOLD([ac_barriers_name_string], [Barriers ${PACKAGE_VERSION}])
AC_STRING_APPEND_FORMAT_BOLD([ac_barriers_final_msg], [You can run 'make', 'make check', and 'make install' now!])

AC_MSG_NOTICE([

======================================
    $ac_barriers_name_string
======================================

Configured successful with the following options:

Extra Features
---------------
  * SECIS Extension           : ${result_secis}

Install Directories
-------------------
  * Executables               : $result__bindir
  * Extra Data                : $result__datadir
  * Man pages                 : $result__mandir

$ac_barriers_final_msg])
