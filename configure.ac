AC_INIT([sisl], 0.1.1)
AC_CONFIG_SRCDIR([src/sisl.h])

# Making releases:
#   SISL_MICRO_VERSION += 1;
#   SISL_INTERFACE_AGE += 1;
#   SISL_BINARY_AGE += 1;
# if any functions have been added, set SISL_INTERFACE_AGE to 0.
# if backwards compatibility has been broken,
# set SISL_BINARY_AGE and SISL_INTERFACE_AGE to 0.
#
SISL_MAJOR_VERSION=0
SISL_MINOR_VERSION=1
SISL_MICRO_VERSION=1
SISL_INTERFACE_AGE=1
SISL_BINARY_AGE=1
SISL_VERSION=$SISL_MAJOR_VERSION.$SISL_MINOR_VERSION.$SISL_MICRO_VERSION

AC_SUBST(SISL_MAJOR_VERSION)
AC_SUBST(SISL_MINOR_VERSION)
AC_SUBST(SISL_MICRO_VERSION)
AC_SUBST(SISL_VERSION)

# libtool versioning
LT_RELEASE=$SISL_MAJOR_VERSION.$SISL_MINOR_VERSION
LT_CURRENT=`expr $SISL_MICRO_VERSION - $SISL_INTERFACE_AGE`
LT_REVISION=$SISL_INTERFACE_AGE
LT_AGE=`expr $SISL_BINARY_AGE - $SISL_INTERFACE_AGE`
AC_SUBST(LT_RELEASE)
AC_SUBST(LT_CURRENT)
AC_SUBST(LT_REVISION)
AC_SUBST(LT_AGE)

# For automake.
VERSION=$SISL_VERSION
PACKAGE=sisl

AM_INIT_AUTOMAKE(subdir-objects)

# Specify a configuration file
AM_CONFIG_HEADER(config.h)

AC_DEFINE_UNQUOTED(SISL_MAJOR_VERSION, $SISL_MAJOR_VERSION, [Major version])
AC_DEFINE_UNQUOTED(SISL_MINOR_VERSION, $SISL_MINOR_VERSION, [Minor version])
AC_DEFINE_UNQUOTED(SISL_MICRO_VERSION, $SISL_MICRO_VERSION, [Micro version])
AC_DEFINE_UNQUOTED(SISL_INTERFACE_AGE, $SISL_INTERFACE_AGE, [Interface age])
AC_DEFINE_UNQUOTED(SISL_BINARY_AGE, $SISL_BINARY_AGE, [Binary age])

dnl Initialize libtool
AM_PROG_LIBTOOL

dnl Initialize maintainer mode
AM_MAINTAINER_MODE

AC_PROG_CC

if test x$GCC = xyes ; then
  CFLAGS="$CFLAGS -Wall -Werror-implicit-function-declaration -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations"
fi

AC_PROG_AWK
AC_ISC_POSIX
AC_HEADER_STDC

AC_CHECK_LIB(m, cos)
AC_CONFIG_MACRO_DIR([m4])

AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_F77_DUMMY_MAIN
AC_F77_WRAPPERS
ACX_BLAS
ACX_LAPACK

blas_libs="$BLAS_LIBS"
lapack_libs="$LAPACK_LIBS"

PKG_CHECK_MODULES([GSL], [gsl]) 

AC_CHECK_PROG(WMPI, wmpi-config, yes, no)
if test x$WMPI = xyes ; then
   AC_DEFINE(HAVE_WMPI, 1, [MPI wrapper library available])
   wmpi_cflags=`wmpi-config --cflags`
   wmpi_libs=`wmpi-config --libs`
   wmpi-config --cc --mpi >& /dev/null
   status=$?
   if test $status -eq 0  ; then
  echo \#################################################################
  echo \#
  echo \# WMPI parallel support: parallel tests will be conducted
  echo \#
  echo \#################################################################
	mpi_tests=mpitests
   else
  echo \#################################################################
  echo \#
  echo \# WMPI has no parallel support: only serial tests will be conducted
  echo \#
  echo \#################################################################
   fi
else
  AC_DEFINE(HAVE_WMPI, 0, [MPI wrapper library available])
  echo \#################################################################
  echo \#
  echo \# WMPI not found: no parallel support is available
  echo \#
  echo \#################################################################
fi

AC_SUBST(gsl_cflags)
AC_SUBST(gsl_libs)
AC_SUBST(wmpi_cflags)
AC_SUBST(wmpi_libs)
AC_SUBST(blas_libs)
AC_SUBST(lapack_libs)
AC_SUBST(mpi_tests)

CFLAGS="$CFLAGS $GSL_CFLAGS $wmpi_cflags"
LDFLAGS="$LDFLAGS $GSL_LIBS"

AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(LDFLAGS)

LIBS="$GSL_LIBS $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

dnl generate `src/sislconfig.h'
AC_CONFIG_COMMANDS([default-1],[[
	echo creating src/sislconfig.h
	cat >src/sislconfig.h <<\__EOF
/* sislconfig.h
 * 
 * This is a generated file.  Please modify `configure.ac'
 */

#ifndef SISLCONFIG_H
#define SISLCONFIG_H

__EOF
	case "$build" in
	*-cray-unicos*)
	  echo "/* CRAY inlining directives */" >> src/sislconfig.h
	  awk '
BEGIN {
  nc = 10;
  nf = 0;
}
{
  if ($1 == "G_INLINE_FUNC") {
    nc = 0;
    start = 2;
  }
  else
    start = 1;
  for (i = start; i <= NF; i++) {
    if (nc == 1) {
      if (substr ($i, 1, 1) != "*") {
	in_list = 0;
	for (j = 0; j < nf && !inlist; j++)
	  if ($i == f[j])
	    inlist = 1;
	if (!inlist)
          f[nf++] = $i;
	nc = 10;
      }
    }
    else
      nc++;
  }
}
END {
  if (nf > 0) {
    printf ("#pragma _CRI inline %s", f[0]);
    for (i = 1; i < nf; i++)
      printf (", %s", f[i]);
    printf ("\n");
  }
}
' < src/sisl.h >> src/sislconfig.h
	;;
	esac
	cat >>src/sislconfig.h <<__EOF
$sisl_defines 

#endif /* SISLCONFIG_H */
__EOF
]],[[
	sisl_defines="
#define SISL_MAJOR_VERSION $SISL_MAJOR_VERSION
#define SISL_MINOR_VERSION $SISL_MINOR_VERSION
#define SISL_MICRO_VERSION $SISL_MICRO_VERSION"
	build="$build"
]])

AC_CONFIG_FILES([
Makefile
src/Makefile
src/sisl-config
tests/Makefile
mpitests/Makefile
mpitests/sisl-test.sh
doc/Makefile
doc/html/Makefile
])

AC_OUTPUT
