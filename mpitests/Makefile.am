## Process this file with automake to produce Makefile.in

CC = `wmpi-config --cc --mpi`
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)
LDADD = $(top_builddir)/src/.libs/libsisl.a 
#AM_CFLAGS = @wmpi_cflags@
#AM_LDFLAGS = @wmpi_libs@
AM_CFLAGS = `wmpi-config --cflags --mpi`
AM_LDFLAGS = `wmpi-config --libs --mpi`
DEPS = $(top_builddir)/src/libsisl.la

bin_PROGRAMS = \
	sisl-test.mpi

sisl_test_mpi_SOURCES = ../tests/sisl-test.c

