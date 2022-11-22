SUBDIRS = fastjet .

if MONOLITHIC_PLUGINS
noinst_LTLIBRARIES = libATLASConePlugin.la
else
lib_LTLIBRARIES = libATLASConePlugin.la
endif
libATLASConePlugin_la_CXXFLAGS = $(AM_CXXFLAGS) -I$(srcdir) -I$(srcdir)/../../include
libATLASConePlugin_la_SOURCES = Jet.cc JetConeFinderTool.cc JetSplitMergeTool.cc ATLASConePlugin.cc CommonUtils.hh JetConeFinderTool.hh JetDistances.hh Jet.hh JetSplitMergeTool.hh LorentzVector.hh
libATLASConePlugin_la_LIBADD = ../../src/libfastjet.la

EXTRA_DIST = makefile.static
