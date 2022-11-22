# shared library

SUBDIRS=fastjet
lib_LTLIBRARIES = libfastjettools.la
libfastjettools_la_CXXFLAGS = $(AM_CXXFLAGS) -I$(srcdir)/../include -I$(srcdir)
libfastjettools_la_LIBADD = ../src/libfastjet.la
libfastjettools_la_SOURCES = \
  BackgroundEstimatorBase.cc \
  JetMedianBackgroundEstimator.cc \
  GridMedianBackgroundEstimator.cc \
  Recluster.cc \
  Filter.cc Pruner.cc Subtractor.cc\
  MassDropTagger.cc RestFrameNSubjettinessTagger.cc CASubJetTagger.cc \
  JHTopTagger.cc TopTaggerBase.cc

# EXTRA_DIST = makefile.static 
