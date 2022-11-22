SWIGINTERFACES = fastjet.i fastjet-doc.i
SWIGDEPENDENCIES = ../include/fastjet/config_auto.h\
    $(top_srcdir)/include/fastjet/config.h\
    $(top_srcdir)/include/fastjet/internal/base.hh\
    $(top_srcdir)/include/fastjet/internal/numconsts.hh\
    $(top_srcdir)/include/fastjet/internal/IsBase.hh\
    $(top_srcdir)/include/fastjet/internal/deprecated.hh\
    $(top_srcdir)/include/fastjet/internal/BasicRandom.hh\
    $(top_srcdir)/include/fastjet/SharedPtr.hh\
    $(top_srcdir)/include/fastjet/LimitedWarning.hh\
    $(top_srcdir)/include/fastjet/Error.hh\
    $(top_srcdir)/include/fastjet/PseudoJetStructureBase.hh\
    $(top_srcdir)/include/fastjet/PseudoJet.hh\
    $(top_srcdir)/include/fastjet/FunctionOfPseudoJet.hh\
    $(top_srcdir)/include/fastjet/RangeDefinition.hh\
    $(top_srcdir)/include/fastjet/Selector.hh\
    $(top_srcdir)/include/fastjet/JetDefinition.hh\
    $(top_srcdir)/include/fastjet/CompositeJetStructure.hh\
    $(top_srcdir)/include/fastjet/WrappedStructure.hh\
    $(top_srcdir)/include/fastjet/ClusterSequenceStructure.hh\
    $(top_srcdir)/include/fastjet/ClusterSequence.hh\
    $(top_srcdir)/include/fastjet/RectangularGrid.hh\
    $(top_srcdir)/include/fastjet/NNBase.hh\
    $(top_srcdir)/include/fastjet/NNH.hh\
    $(top_srcdir)/include/fastjet/NNFJN2Plain.hh\
    $(top_srcdir)/include/fastjet/NNFJN2Tiled.hh\
    $(top_srcdir)/include/fastjet/GhostedAreaSpec.hh\
    $(top_srcdir)/include/fastjet/AreaDefinition.hh\
    $(top_srcdir)/include/fastjet/ClusterSequenceAreaBase.hh\
    $(top_srcdir)/include/fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh\
    $(top_srcdir)/include/fastjet/ClusterSequenceActiveArea.hh\
    $(top_srcdir)/include/fastjet/ClusterSequence1GhostPassiveArea.hh\
    $(top_srcdir)/include/fastjet/ClusterSequencePassiveArea.hh\
    $(top_srcdir)/include/fastjet/ClusterSequenceVoronoiArea.hh\
    $(top_srcdir)/include/fastjet/ClusterSequenceArea.hh\
    $(top_srcdir)/tools/fastjet/tools/Transformer.hh\
    $(top_srcdir)/tools/fastjet/tools/Boost.hh\
    $(top_srcdir)/tools/fastjet/tools/Recluster.hh\
    $(top_srcdir)/tools/fastjet/tools/Filter.hh\
    $(top_srcdir)/tools/fastjet/tools/Pruner.hh\
    $(top_srcdir)/tools/fastjet/tools/CASubJetTagger.hh\
    $(top_srcdir)/tools/fastjet/tools/MassDropTagger.hh\
    $(top_srcdir)/tools/fastjet/tools/RestFrameNSubjettinessTagger.hh\
    $(top_srcdir)/tools/fastjet/tools/TopTaggerBase.hh\
    $(top_srcdir)/tools/fastjet/tools/JHTopTagger.hh\
    $(top_srcdir)/tools/fastjet/tools/BackgroundEstimatorBase.hh\
    $(top_srcdir)/tools/fastjet/tools/JetMedianBackgroundEstimator.hh\
    $(top_srcdir)/tools/fastjet/tools/GridMedianBackgroundEstimator.hh\
    $(top_srcdir)/tools/fastjet/tools/Subtractor.hh\
    $(srcdir)/FastJetPythonExtensions.hh

# include the interface file for completeness
EXTRA_DIST=$(SWIGINTERFACES) FastJetPythonExtensions.hh

if BUILD_INCLUDES_SWIG
BUILT_SOURCES = swig_wrap.cpp fastjet.py
swig_wrap.cpp : $(SWIGINTERFACES) $(SWIGDEPENDENCIES)
	$(SWIG) $(AX_SWIG_PYTHON_OPT) -I../include -I$(top_srcdir)/include -I$(top_srcdir)/tools -o $@ $<
fastjet.py: swig_wrap.cpp
	@echo "fastjet.py is built at the same time as swig_wrap.cpp"
else
# you're responsible for making sure this is available (through the
# distrib or by running swig for developpers)
EXTRA_DIST+=swig_wrap.cpp fastjet.py
endif

python_PYTHON = fastjet.py
pyexec_LTLIBRARIES = _fastjet.la
_fastjet_la_LIBTOOLFLAGS = --silent
_fastjet_la_SOURCES = swig_wrap.cpp
_fastjet_la_CPPFLAGS = $(PYTHON_INCLUDE) $(AM_CXXFLAGS) $(CGAL_CPPFLAGS) -I../include -I$(top_srcdir)/include -I$(top_srcdir)/tools -DSWIG -DSWIG_TYPE_TABLE=fastjet
_fastjet_la_LDFLAGS = -module
_fastjet_la_LIBADD = ../src/libfastjet.la ../tools/libfastjettools.la $(CGAL_LIBS)

# include the interface file for completeness
pyinterfacedir=$(datadir)/fastjet/pyinterface
pyinterface_DATA=$(SWIGINTERFACES) FastJetPythonExtensions.hh

# this gets evaluated before install-exec
# we use that to silence litool
install-exec-local:
	$(eval AM_LIBTOOLFLAGS = "--silent")

# this gets evaluated after install-exec
# we use that restore the flags to nothing (as they're supposed to be)
install-exec-hook:
	$(eval AM_LIBTOOLFLAGS = "")
	@echo "------------------------------------------------------------"
	@echo "The fastjet.py python module is installed in "
	@echo "  $(DESTDIR)$(pythondir)"
	@echo "and the _fastjet.so library is installed in" 
	@echo "  $(DESTDIR)$(pyexecdir)"
	@echo "You should make sure that these are in your PYTHONPATH"
	@echo "------------------------------------------------------------"


#----------------------------------------------------------------------
# If we were to let the end-user run swig, then we should have the following:

#BUILT_SOURCES = $(srcdir)/swig_wrap.cpp
#SWIG_SOURCES = fastjet.i
#
#pkgpython_PYTHON = fastjet.py
#pkgpyexec_LTLIBRARIES = _fastjet.la
#_fastjet_la_SOURCES = swig_wrap.cpp
#_fastjet_la_CPPFLAGS = $(AX_SWIG_PYTHON_CPPFLAGS) $(PYTHON_INCLUDE) $(AM_CXXFLAGS) $(CGAL_CPPFLAGS) -I../include -I$(top_srcdir)/include -I$(top_srcdir)/tools -DSWIG
#_fastjet_la_LDFLAGS = -module
#_fastjet_la_LIBADD = ../src/libfastjet.la ../tools/libfastjettools.la
#
#$(srcdir)/swig_wrap.cpp : $(SWIG_SOURCES)
#	$(SWIG) $(AX_SWIG_PYTHON_OPT) -I$(top_srcdir)/include -I$(top_srcdir)/tools -o $@ $<

