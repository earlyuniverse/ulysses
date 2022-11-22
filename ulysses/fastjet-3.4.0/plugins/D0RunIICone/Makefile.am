SUBDIRS = fastjet

if MONOLITHIC_PLUGINS
noinst_LTLIBRARIES = libD0RunIIConePlugin.la
else
lib_LTLIBRARIES = libD0RunIIConePlugin.la
endif
libD0RunIIConePlugin_la_CXXFLAGS = $(D0_CXXFLAGS) -I$(srcdir) -I$(srcdir)/../../include
libD0RunIIConePlugin_la_SOURCES = D0RunIIConePlugin.cc
libD0RunIIConePlugin_la_LIBADD = ../../src/libfastjet.la

#EXTRA_DIST = makefile.static pxcone.h
EXTRA_DIST = ConeJetInfo.hpp  ConeSplitMerge.hpp  ILConeAlgorithm.hpp  ProtoJet.hpp HepEntity.h  inline_maths.h

