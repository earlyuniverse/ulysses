SUBDIRS = yoda

EXTRA_DIST = yoda1/__init__.py

if ENABLE_PYEXT
PYEXT_ENV = CC="$(CXX)" CXX="$(CXX)" CXXFLAGS="$(PYEXT_CXXFLAGS)" ARCHFLAGS=""

if ENABLE_ROOT
  PYEXT_ENV += BUILD_ROOTCOMPAT="yes"
endif

## Always force setup.py, it's not good at guessing what needs to rebuild
all-local:
	cd $(srcdir); \
	$(PYEXT_ENV) YODA_LOCAL=1 $(PYTHON) \
	  $(abs_builddir)/setup.py build \
        --build-base $(abs_builddir)/build \
        --verbose \
        --force

install-exec-local:
	$(PYEXT_ENV) $(PYTHON) \
	  $(abs_builddir)/setup.py install \
	    --root "$(DESTDIR)///" \
	    --prefix "$(prefix)" \
	    --verbose \
	    --skip-build \
	    --force 

uninstall-local:
	dir1=$(DESTDIR)$(pythondir); \
	dir2=$${dir1/lib/lib64}; \
	for i in $$dir1 $$dir2; do \
	 rm -rf $$i/yoda; \
	 rm -rf $$i/yoda1; \
	 rm -rf $$i/yoda*.egg-info; \
	done

clean-local:
	rm -rf $(builddir)/build

endif
