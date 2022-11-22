SUBDIRS = fastjet .

if MONOLITHIC_PLUGINS
noinst_LTLIBRARIES = libTrackJetPlugin.la
else
lib_LTLIBRARIES = libTrackJetPlugin.la
endif
libTrackJetPlugin_la_CXXFLAGS = $(AM_CXXFLAGS) -I$(srcdir) -I$(srcdir)/../../include
libTrackJetPlugin_la_SOURCES = TrackJetPlugin.cc
libTrackJetPlugin_la_LIBADD = ../../src/libfastjet.la

EXTRA_DIST = makefile.static
