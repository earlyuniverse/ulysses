SUBDIRS = fastjet .

if MONOLITHIC_PLUGINS
noinst_LTLIBRARIES = libGridJetPlugin.la
else
lib_LTLIBRARIES = libGridJetPlugin.la
endif
libGridJetPlugin_la_CXXFLAGS = ${AM_CXXFLAGS} -I$(srcdir) -I$(srcdir)/../../include
libGridJetPlugin_la_SOURCES = GridJetPlugin.cc
libGridJetPlugin_la_LIBADD = ../../src/libfastjet.la

EXTRA_DIST = makefile.static
