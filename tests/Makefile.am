## Process this file with automake to produce Makefile.in

CC = `wmpi-config --cc`
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)
LDADD = $(top_builddir)/src/.libs/libsisl.a 

AM_CFLAGS = @wmpi_cflags@
AM_LDFLAGS = @wmpi_libs@

DEPS = $(top_builddir)/src/libsisl.la

bin_PROGRAMS = \
	sisl-test
