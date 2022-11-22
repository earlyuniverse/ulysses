# we only install the 'siscone' and 'area' programs
# since the remaining ones are samples and test programs
# we leave them in the 'examples' folder.
noinst_PROGRAMS = siscone siscone_area test times sample spherical

EXTRA_DIST = options.h makefile.static

SUBDIRS = events

# specify sources for each program to build
siscone_SOURCES = options.cpp main.cpp
siscone_CXXFLAGS = $(AM_CXXFLAGS) -I$(srcdir)/.. -I..
siscone_LDADD = ../siscone/libsiscone.la

siscone_area_SOURCES = options.cpp area.cpp
siscone_area_CXXFLAGS = $(AM_CXXFLAGS) -I$(srcdir)/.. -I..
siscone_area_LDADD = ../siscone/libsiscone.la

test_SOURCES = test.cpp
test_CXXFLAGS = $(AM_CXXFLAGS) -I$(srcdir)/.. -I..
test_LDADD = ../siscone/libsiscone.la

times_SOURCES = times.cpp
times_CXXFLAGS = $(AM_CXXFLAGS) -I$(srcdir)/.. -I..
times_LDADD = ../siscone/libsiscone.la

sample_SOURCES = sample.cpp
sample_CXXFLAGS = $(AM_CXXFLAGS) -I$(srcdir)/.. -I..
sample_LDADD = ../siscone/libsiscone.la

spherical_SOURCES = spherical.cpp
spherical_CXXFLAGS = $(AM_CXXFLAGS) -I$(srcdir)/.. -I..
spherical_LDADD = ../siscone/spherical/libsiscone_spherical.la\
	../siscone/libsiscone.la

