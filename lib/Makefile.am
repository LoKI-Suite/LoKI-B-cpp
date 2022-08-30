SUBDIRS = .

include $(top_srcdir)/Makefile.rules

VPATH=$(top_srcdir)/source

libloki_b_la_LDFLAGS = -no-undefined
libloki_b_la_SOURCES = \
	CrossSection.cpp \
	EedfCollisions.cpp \
	ElectronKinetics.cpp \
	EedfMixture.cpp \
	Enumeration.cpp \
	Gas.cpp \
	GasMixture.cpp \
	Gnuplot.cpp \
	Grid.cpp \
	JobSystem.cpp \
	LinearAlgebra.cpp \
	LookupTable.cpp \
	Output.cpp \
	Parse.cpp \
	Power.cpp \
	Setup.cpp \
	Simulation.cpp \
	StateEntry.cpp \
	WorkingConditions.cpp \
	json.cpp

libdir=${prefix}/lib
lib_LTLIBRARIES = \
	libloki-b.la
