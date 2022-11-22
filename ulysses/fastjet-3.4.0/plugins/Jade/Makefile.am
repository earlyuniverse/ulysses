SUBDIRS = fastjet .

if MONOLITHIC_PLUGINS
noinst_LTLIBRARIES = libJadePlugin.la
else
lib_LTLIBRARIES = libJadePlugin.la
endif
libJadePlugin_la_CXXFLAGS = $(AM_CXXFLAGS) -I$(srcdir) -I$(srcdir)/../../include
libJadePlugin_la_SOURCES = JadePlugin.cc
libJadePlugin_la_LIBADD = ../../src/libfastjet.la

EXTRA_DIST = makefile.static
