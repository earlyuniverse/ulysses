SUBDIRS = tinyxml yamlcpp .

lib_LTLIBRARIES = libYODA.la

libYODA_la_SOURCES = \
    Exceptions.cc \
    Reader.cc \
    ReaderYODA.cc \
    ReaderFLAT.cc \
    ReaderAIDA.cc \
    Writer.cc \
    WriterYODA.cc \
    WriterFLAT.cc \
    WriterAIDA.cc \
    Dbn0D.cc \
    Dbn1D.cc \
    Counter.cc \
    Histo1D.cc \
    Histo2D.cc \
    Profile1D.cc \
    Profile2D.cc \
    Scatter1D.cc \
    Scatter2D.cc \
    Scatter3D.cc \
    Point1D.cc \
    Point2D.cc \
    Point3D.cc

libYODA_la_LDFLAGS = -avoid-version
libYODA_la_LIBADD = $(builddir)/tinyxml/libyoda-tinyxml.la $(builddir)/yamlcpp/libyoda-yaml-cpp.la
libYODA_la_CPPFLAGS = $(AM_CPPFLAGS) -DTIXML_USE_STL -I$(srcdir)/yamlcpp -I$(srcdir) -DYAML_NAMESPACE=YODA_YAML

EXTRA_DIST = zstr
