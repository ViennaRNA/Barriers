#!/bin/sh
# 
# Run this before configure
#
# This file blatantly ripped off from subversion.
#
# $Id: autogen.sh,v 1.3 2008/01/09 10:14:00 xtof Exp $

set -e
# Produce aclocal.m4, so autoconf gets the automake macros it needs
echo "Creating aclocal.m4..."
aclocal

autoheader

# Produce all the `Makefile.in's, verbosely, and create neat missing things
# like `libtool', `install-sh', etc.
automake --add-missing --gnu --copy

# If there's a config.cache file, we may need to delete it.  
# If we have an existing configure script, save a copy for comparison.
if [ -f config.cache ] && [ -f configure ]; then
  cp configure configure.$$.tmp
fi

# Produce ./configure
echo "Creating configure..."
autoconf

echo ""
echo "You can run ./configure now."
echo ""
