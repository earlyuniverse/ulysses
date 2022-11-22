SUBDIRS = fastjet .

if MONOLITHIC_PLUGINS
noinst_LTLIBRARIES = libEECambridgePlugin.la
else
lib_LTLIBRARIES = libEECambridgePlugin.la
endif
libEECambridgePlugin_la_CXXFLAGS = $(AM_CXXFLAGS) -I$(srcdir) -I$(srcdir)/../../include
libEECambridgePlugin_la_SOURCES = EECambridgePlugin.cc
libEECambridgePlugin_la_LIBADD = ../../src/libfastjet.la

EXTRA_DIST = makefile.static
