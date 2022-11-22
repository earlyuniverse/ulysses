# build information for the SISCone library
# this is built as a libtool lib.
lib_LTLIBRARIES = libsiscone_spherical.la

# On top of the configured flags, we need to access the the 
# top dir for includes of the form <siscone/*.h> to work
# Note that we use "siscone/..." for clarity.
# Then we need to access the config.h file which is one step away
# since we use <siscone/config.h> (for consistency), this means ../..
libsiscone_spherical_la_CXXFLAGS = $(AM_CXXFLAGS) -I$(srcdir)/../.. -I../..
libsiscone_spherical_la_SOURCES = geom_2d.cpp momentum.cpp hash.cpp\
   vicinity.cpp protocones.cpp split_merge.cpp siscone.cpp
libsiscone_spherical_la_LIBADD = ../libsiscone.la

EXTRA_DIST = makefile.static

# install the SISCone headers
sisconeincludedir = $(includedir)/siscone/spherical
sisconeinclude_HEADERS = geom_2d.h\
  hash.h\
  momentum.h\
  protocones.h\
  siscone.h\
  split_merge.h\
  vicinity.h
