SUBDIRS = siscone fastjet .

if MONOLITHIC_PLUGINS
noinst_LTLIBRARIES = libSISConePlugin.la
else
lib_LTLIBRARIES = libSISConePlugin.la
endif
libSISConePlugin_la_CXXFLAGS = $(AM_CXXFLAGS) -I$(srcdir) -I$(srcdir)/siscone -I$(srcdir)/../../include -I./siscone
libSISConePlugin_la_SOURCES = SISConeBasePlugin.cc SISConePlugin.cc SISConeSphericalPlugin.cc
libSISConePlugin_la_LIBADD = ../../src/libfastjet.la siscone/siscone/libsiscone.la siscone/siscone/spherical/libsiscone_spherical.la

EXTRA_DIST = makefile.static
