SUBDIRS = .

include $(top_srcdir)/Makefile.rules

VPATH=$(top_srcdir)/loki-b/source

libloki_b_la_SOURCES = \
	CrossSection.cpp \
	EedfCollision.cpp \
	EedfGas.cpp \
	EedfGasMixture.cpp \
	ElectronKinetics.cpp \
	Enumeration.cpp \
	GasBase.cpp \
	GasMixtureBase.cpp \
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
