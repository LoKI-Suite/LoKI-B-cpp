ACLOCAL_AMFLAGS=-I autodir

SUBDIRS = lib . testsuite @DOC_SUBDIR@

include $(top_srcdir)/Makefile.rules

VPATH=$(top_srcdir)/loki-b/source

loki_SOURCES = main.cpp
loki_LDADD = $(top_builddir)/lib/libloki-b.la

bin_PROGRAMS = loki
