SUBDIRS = siscone examples
EXTRA_DIST = makefile.static Doxyfile

# local macros for aclocal
ACLOCAL_AMFLAGS = -I m4

# AX_PREFIX_CONFIG_H leaves some files undeleted after a make distclean
#
# we could use DISTCLEAN=... here but this would also overwrite the defaults
# So, instead, we use the (apparently recommended) -local target to extend 
# the default behaviour
distclean-local: distclean-ax-prefix-config-h
distclean-ax-prefix-config-h:
	rm -f _configs.sed siscone/config.h
