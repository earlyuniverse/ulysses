SUBDIRS = CDFcode fastjet .

if MONOLITHIC_PLUGINS
noinst_LTLIBRARIES = libCDFConesPlugin.la
else
lib_LTLIBRARIES = libCDFConesPlugin.la
endif
libCDFConesPlugin_la_CXXFLAGS = $(AM_CXXFLAGS) -I$(srcdir) -I$(srcdir)/../../include -I$(srcdir)/CDFcode
libCDFConesPlugin_la_SOURCES = CDFJetCluPlugin.cc CDFMidPointPlugin.cc
libCDFConesPlugin_la_LIBADD = CDFcode/libCDFcode.la ../../src/libfastjet.la


EXTRA_DIST = makefile.static
