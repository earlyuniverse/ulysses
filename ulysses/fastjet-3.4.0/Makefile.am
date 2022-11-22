# redirect to subdirectories
SUBDIRS = src include plugins tools example fortran_wrapper doc
if BUILD_PYEXT
SUBDIRS += pyinterface
endif

# install the fastjet-config script in the bin directory
bin_SCRIPTS = fastjet-config

# tests done using make check
TEST_ENVIRONMENT = $(SHELL)
TESTS = test-compare.sh

# local macros for aclocal
ACLOCAL_AMFLAGS = -I m4

# additional files we need in the distribution
EXTRA_DIST = test-compare.sh makefile.static test-static.sh test-script-output-orig.txt Doxyfile BUGS

# AX_PREFIX_CONFIG_H leaves some files undeleted after a make distclean
#
# we could use DISTCLEAN=... here but this would also overwrite the defaults
# So, instead, we use the (apparently recommended) -local target to extend 
# the default behaviour
distclean-local: distclean-ax-prefix-config-h
distclean-ax-prefix-config-h:
	rm -f _configs.sed include/fastjet/config_auto.h
