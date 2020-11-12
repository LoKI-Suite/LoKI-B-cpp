ACLOCAL_AMFLAGS=-I autodir

SUBDIRS = lib .

AM_CPPFLAGS = -I$(top_srcdir)/loki-b -I$(top_srcdir)/include
AM_CXXFLAGS = -g -pedantic -Wall -std=c++17

#all: loki

VPATH=$(top_srcdir)/loki-b/source

loki_SOURCES = main.cpp
loki_LDADD = $(top_builddir)/lib/libloki-b.la

bin_PROGRAMS = loki
