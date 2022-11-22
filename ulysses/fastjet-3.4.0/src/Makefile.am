# shared library

lib_LTLIBRARIES = libfastjet.la
libfastjet_la_CXXFLAGS = $(AM_CXXFLAGS) -I$(srcdir)/../include
libfastjet_la_SOURCES = DnnPlane.cc\
	Dnn4piCylinder.cc Dnn3piCylinder.cc Dnn2piCylinder.cc\
	ClusterSequence.cc PseudoJet.cc Selector.cc\
	ClusterSequence_N2.cc ClusterSequence_TiledN2.cc\
	ClusterSequence_Delaunay.cc ClusterSequence_DumbN3.cc\
	ClusterSequence_CP2DChan.cc ClosestPair2D.cc MinHeap.cc\
	ClusterSequenceAreaBase.cc ClusterSequenceActiveAreaExplicitGhosts.cc\
	ClusterSequenceArea.cc GhostedAreaSpec.cc ClusterSequenceActiveArea.cc\
	Voronoi.cc ClusterSequenceVoronoiArea.cc\
	ClusterSequencePassiveArea.cc ClusterSequence1GhostPassiveArea.cc\
	PseudoJetStructureBase.cc ClusterSequenceStructure.cc\
	BasicRandom.cc JetDefinition.cc Error.cc AreaDefinition.cc\
	RangeDefinition.cc CompositeJetStructure.cc \
	FunctionOfPseudoJet.cc\
	LimitedWarning.cc\
	TilingExtent.cc \
	LazyTiling9.cc \
        LazyTiling9Alt.cc \
        LazyTiling25.cc \
        LazyTiling9SeparateGhosts.cc \
        RectangularGrid.cc

EXTRA_DIST = makefile.static genconfig.sh
