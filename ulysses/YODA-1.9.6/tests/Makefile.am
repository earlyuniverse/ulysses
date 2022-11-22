TEST_EXTENSIONS = .py .sh

PYTESTS = \
  test-unit.py \
  test-counter.py \
  test-h1d.py \
  test-p1d.py \
  test-h2d.py \
  test-p2d.py \
  test-s1d.py \
  test-s2d.py \
  test-div.py \
  test-rebin.py \
  test-io.py \
  test-iofilter.py \
  test-operators.py

SHTESTS = \
  test-yodamerge.sh \
  test-yodahist.sh \
  test-yodals.sh \
  test-yodacmp.sh \
  test-yoda2flat.sh \
  test-yodascale.sh \
  test-yodaplot.sh \
  test-yodacnv.sh \
  test-yoda2yoda.sh

EXTRA_DIST = $(PYTESTS) $(SHTESTS) \
  testreader.sh testwriter.sh \
  test-yoda2root.sh \
  test.yoda test.yoda.gz \
  test1.yoda test2.yoda \
  iofilter.yoda \
  rivetexample.yoda \
  yodals-ref.txt \
  yodahist-ref.yoda \
  yodahist-fills.txt \
  yoda2flat-ref.dat \
  merged12-ref.yoda merged12pi-ref.yoda

check_PROGRAMS = \
  testtraits \
  testannotations \
  testweights \
  testbinsearcher \
  testwriter \
  testreader \
  testhisto1Da testhisto1Db \
  testhisto2Da \
  testprofile1Da \
  testindexedset \
  testsortedvector\
  testhisto1Dcreate \
  testhisto1Dfill \
  testhisto1Dmodify \
  testprofile1Dcreate \
  testprofile1Dfill \
  testprofile1Dmodify \
  testscatter2Dcreate \
  testscatter2Dmodify \
  testhisto2Dcreate

#  testhisto2Dfill \
#  testhisto2Dmodify \
#  testscatter3Dcreate \
#  testscatter3Dmodify

AM_LDFLAGS = -L$(top_builddir)/src -lYODA
AM_CPPFLAGS += -DYODA_TESTS_SRC="\"$(srcdir)\""

testtraits_SOURCES = TestTraits.cc
testannotations_SOURCES = TestAnnotations.cc
testweights_SOURCES = TestWeights.cc
testbinsearcher_SOURCES = TestBinSearcher.cc
testwriter_SOURCES = TestWriter.cc
testreader_SOURCES = TestReader.cc
testhisto1Da_SOURCES = TestHisto1Da.cc
testhisto1Db_SOURCES = TestHisto1Db.cc
testprofile1Da_SOURCES = TestProfile1Da.cc
testindexedset_SOURCES = TestIndexedSet.cc
testsortedvector_SOURCES = TestSortedVector.cc
testhisto1Dcreate_SOURCES = Histo1D/H1DCreate.cc
testhisto1Dfill_SOURCES = Histo1D/H1DFill.cc
testhisto1Dmodify_SOURCES = Histo1D/H1DModify.cc
testprofile1Dcreate_SOURCES = Profile1D/P1DCreate.cc
testprofile1Dfill_SOURCES = Profile1D/P1DFill.cc
testprofile1Dmodify_SOURCES = Profile1D/P1DModify.cc
testhisto2Da_SOURCES = TestHisto2Da.cc
# testhisto2De_SOURCES = TestHisto2Derase.cc
testhisto2Dcreate_SOURCES = Histo2D/H2DCreate.cc
#testhisto2Dfill_SOURCES = Histo2D/H2DFill.cc
#testhisto2Dmodify_SOURCES = Histo2D/H2DModify.cc
testscatter2Dcreate_SOURCES = Scatter2D/S2DCreate.cc
testscatter2Dmodify_SOURCES = Scatter2D/S2DModify.cc
# testpoint3D_SOURCES = TestPoint3D.cc
# testscatter3D_SOURCES = TestScatter3D.cc
# testscatter3Dcreate_SOURCES = Scatter3D/S3DCreate.cc
# testscatter3Dmodify_SOURCES = Scatter3D/S3DModify.cc


TESTS_ENVIRONMENT = \
  LD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(LD_LIBRARY_PATH) \
  DYLD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(DYLD_LIBRARY_PATH) \
  PYTHONPATH=$(top_builddir)/pyext/build/$(YODA_SETUP_PY_PATH):$(PYTHONPATH) \
  PATH=$(top_builddir)/bin:$(top_srcdir)/bin:$(PATH) \
  YODA_TESTS_SRC=$(srcdir)


TESTS = \
  testtraits \
  testannotations \
  testweights \
  testbinsearcher \
  testwriter.sh \
  testreader.sh \
  testhisto1Da \
  testhisto1Db \
  testhisto2Da \
  testprofile1Da \
  testindexedset \
  testsortedvector \
  testhisto1Dcreate \
  testhisto1Dfill \
  testhisto1Dmodify \
  testprofile1Dcreate \
  testprofile1Dfill \
  testprofile1Dmodify \
  testscatter2Dcreate \
  testscatter2Dmodify \
  testhisto2Dcreate

testreader.log: testwriter.log

TESTS += $(PYTESTS) $(SHTESTS)

CLEANFILES = \
  h1d.yoda h1d.dat \
  p1d.yoda p1d.dat \
  h2d.yoda h2d.dat \
  p2d.yoda p2d.dat \
  s1d.yoda s2d.yoda \
  testwriter1.yoda testwriter2.yoda testwriter2.yoda.gz \
  foo_bar_baz.dat \
  counter.yoda \
  test.aida \
  y2y_3.yoda

if ENABLE_ROOT
  TESTS += test-yoda2root.sh
  CLEANFILES += test1.root
endif
