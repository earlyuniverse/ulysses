ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src include tests bin

if ENABLE_PYEXT
SUBDIRS += pyext
endif

EXTRA_FILES = main.dox PLOTKEYS

pkgconfigdir = $(libdir)/pkgconfig
dist_pkgconfig_DATA = yoda.pc

## Deal with the Doxygen stuff
# TODO: Improve, without the Doxy performance penalty on every call of "make"

if WITH_DOXYGEN

.PHONY: dox
dox:
	@echo "Running Doxygen..."
	$(DOXYGEN) $(top_srcdir)/Doxyfile

.PHONY: pydoc
pydoc:
	@echo "Running pydoc..."
	cd pydoc && $(MKDIR_P) _static && $(MAKE) html

all-local: #dox
	@true

#install-data-local:
#	$(install_sh_DATA) doxy/html/* $(DESTDIR)$(pkgdatadir)/doc/html

mostlyclean-local:
	rm -rf $(top_builddir)/doxy

uninstall-local:
	rm -rf $(DESTDIR)$(pkgdatadir)/doc
	test ! -d $(DESTDIR)$(pkgdatadir) || rmdir $(DESTDIR)$(pkgdatadir) || true

endif
