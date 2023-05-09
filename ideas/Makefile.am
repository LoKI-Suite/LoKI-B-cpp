SUBDIRS = .

include $(top_srcdir)/Makefile.rules

json2lxcat_SOURCES = json2lxcat.cpp
json2lxcat_LDADD =

test_statecontainer_SOURCES = test_statecontainer.cpp
test_statecontainer_LDADD = $(LOKIB_LIBS)

bin_PROGRAMS = \
	json2lxcat \
	test_statecontainer

TESTS = \
	test_statecontainer
