AM_CPPFLAGS = -I$(top_srcdir)/loki-b -I$(top_srcdir)/include
AM_CXXFLAGS = -g -pedantic -Wall -std=c++17
AM_LDFLAGS = -no-undefined

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
