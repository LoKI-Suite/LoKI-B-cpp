SUBDIRS = .

include $(top_srcdir)/Makefile.rules

VPATH=$(top_srcdir)/loki-b/tests

bolsig_extract_SOURCES = bolsig_extract.cpp

compare_tables_SOURCES = compare_tables.cpp
compare_tables_LDADD = $(top_builddir)/lib/libloki-b.la

druyvesteyn_SOURCES = druyvesteyn.cpp
druyvesteyn_LDADD = $(top_builddir)/lib/libloki-b.la

test_lut_SOURCES = test_lut.cpp
test_lut_LDADD = $(top_builddir)/lib/libloki-b.la

test_parser_SOURCES = test_parser.cpp
test_parser_LDADD = $(top_builddir)/lib/libloki-b.la

test_reaction_format_SOURCES = test_reaction_format.cpp
test_reaction_format_LDADD = $(top_builddir)/lib/libloki-b.la

test_statebase_SOURCES = test_statebase.cpp
test_statebase_LDADD = $(top_builddir)/lib/libloki-b.la

bin_PROGRAMS = \
	bolsig_extract \
	compare_tables \
	druyvesteyn \
	test_lut \
	test_parser \
	test_reaction_format \
	test_statebase

TESTS = \
	druyvesteyn \
	test_lut \
	test_parser \
	test_reaction_format \
	test_statebase
