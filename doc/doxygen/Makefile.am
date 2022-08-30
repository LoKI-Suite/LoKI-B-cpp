SUBDIRS =

if HAVE_DOXYGEN
doxyfile.stamp:
	$(DOXYGEN) loki-b.doxy
	echo Timestamp > doxyfile.stamp

CLEANFILES = doxyfile.stamp

all-local: doxyfile.stamp

endif 
