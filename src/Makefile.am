## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(includedir) -DG_LOG_DOMAIN=\"SISL\"

bin_SCRIPTS=sisl-config

BUILT_SOURCES= \
	sisl-config \
	sislconfig.h

sisl-config: sisl-config.in

lib_LTLIBRARIES = libsisl.la

libsisl_la_LDFLAGS = -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE) -export-dynamic

libsisl_la_SOURCES = \
	vector.c \
	matrix.c \
	sparse.c \
	iter-solver.c \
	iter-solver-c.c \
	iter-solver.h \
	sisl-logging.c \
	sisl-private.h \
	diagonal.c

include_HEADERS = \
	sisl.h \
	sislconfig.h

CLEANFILES = $(BUILT_SOURCES)
