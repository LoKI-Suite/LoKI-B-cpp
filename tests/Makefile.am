SUBDIRS = .

include $(top_srcdir)/Makefile.rules

bolsig_extract_SOURCES = bolsig_extract.cpp

compare_tables_SOURCES = compare_tables.cpp
compare_tables_LDADD = $(LOKIB_LIBS)

test_druyvesteyn_SOURCES = test_druyvesteyn.cpp
test_druyvesteyn_LDADD = $(LOKIB_LIBS)

test_ee_SOURCES = test_ee.cpp
test_ee_LDADD = $(LOKIB_LIBS)

test_ee_relaxation_SOURCES = test_ee_relaxation.cpp
test_ee_relaxation_LDADD = $(LOKIB_LIBS)

test_events_SOURCES = test_events.cpp
test_events_LDADD = $(LOKIB_LIBS)

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
	test_ee \
	test_ee_relaxation \
	test_events \
	test_lut \
	test_parser \
	test_reaction_format \
	test_statebase

TESTS = \
	test_druyvesteyn \
	test_ee \
	test_events \
	test_lut \
	test_parser \
	test_reaction_format \
	test_statebase
