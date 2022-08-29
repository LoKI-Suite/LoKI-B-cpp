SUBDIRS = .

include $(top_srcdir)/Makefile.rules

bolsig_extract_SOURCES = bolsig_extract.cpp

compare_tables_SOURCES = compare_tables.cpp
compare_tables_LDADD = $(LOKIB_LIBS)

test_druyvesteyn_SOURCES = test_druyvesteyn.cpp
test_druyvesteyn_LDADD = $(LOKIB_LIBS)

test_lut_SOURCES = test_lut.cpp
test_lut_LDADD = $(LOKIB_LIBS)

test_parser_SOURCES = test_parser.cpp
test_parser_LDADD = $(LOKIB_LIBS)

test_reaction_format_SOURCES = test_reaction_format.cpp
test_reaction_format_LDADD = $(LOKIB_LIBS)

test_statebase_SOURCES = test_statebase.cpp
test_statebase_LDADD = $(LOKIB_LIBS)

bin_PROGRAMS = \
	bolsig_extract \
	compare_tables \
	test_druyvesteyn \
	test_lut \
	test_parser \
	test_reaction_format \
	test_statebase

TESTS = \
	test_druyvesteyn \
	test_lut \
	test_parser \
	test_reaction_format \
	test_statebase
