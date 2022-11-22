SUBDIRS = fastjet .

if MONOLITHIC_PLUGINS
noinst_LTLIBRARIES = libNestedDefsPlugin.la
else
lib_LTLIBRARIES = libNestedDefsPlugin.la
endif
libNestedDefsPlugin_la_CXXFLAGS = $(AM_CXXFLAGS) -I$(srcdir) -I$(srcdir)/../../include
libNestedDefsPlugin_la_SOURCES = NestedDefsPlugin.cc
libNestedDefsPlugin_la_LIBADD = ../../src/libfastjet.la

EXTRA_DIST = makefile.static
