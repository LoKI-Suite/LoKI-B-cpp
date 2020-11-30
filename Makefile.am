ACLOCAL_AMFLAGS=-I autodir

SUBDIRS = lib . testsuite @DOC_SUBDIR@

include $(top_srcdir)/Makefile.rules

VPATH=$(top_srcdir)/loki-b/source

loki_SOURCES = main.cpp
loki_LDADD = $(top_builddir)/lib/libloki-b.la

bin_PROGRAMS = loki

CROSSTOOLDIR=$(realpath $(shell dirname $(shell which $(CXX)))/..)
CROSSDLLDIR=$(CROSSTOOLDIR)/$(host_triplet)/bin

# install the MXE libraries. Call this for an MXE stand-alone cross-build.
# it copies the required mxe libs to the bin/ installation subdirectory
install-mxe-libs:
	@echo "Copying MXE's dll libraries from directory $(CROSSDLLDIR) into $(DESTDIR)$(bindir)..."
	@cp -v $(CROSSDLLDIR)/libgcc_s_seh-1.dll $(DESTDIR)$(bindir)
	@cp -v $(CROSSDLLDIR)/libwinpthread-1.dll $(DESTDIR)$(bindir)
	@cp -v $(CROSSDLLDIR)/libstdc++-6.dll $(DESTDIR)$(bindir)
