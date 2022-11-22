.NOTPARALLEL:
	@true

STATIC_SRC = \
    __init__.py \
    axis.pyx \
    errors.cpp errors.hh \
    pyroot_helpers.hh merge.hh \
    rootcompat.pxd rootcompat.pyx \
    search.py plotting.py root.py \
    script_helpers.py

CYTHON_SRC = \
    core.pyx \
    declarations.pxd \
    util.pyx \
    util.pxd \
    include/AnalysisObject.pyx \
    include/Counter.pyx \
    include/Axis1D_BIN1D_DBN.pyx \
    include/Axis2D_BIN2D_DBN.pyx \
    include/Bin.pyx \
    include/Bin1D_DBN.pyx \
    include/Bin2D_DBN.pyx \
    include/Dbn0D.pyx \
    include/Dbn1D.pyx \
    include/Dbn2D.pyx \
    include/Dbn3D.pyx \
    include/Errors.pyx \
    include/Functions.pyx \
    include/Histo1D.pyx \
    include/Histo2D.pyx \
    include/HistoBin1D.pyx \
    include/HistoBin2D.pyx \
    include/IO.pyx \
    include/Point.pyx \
    include/Point1D.pyx \
    include/Point2D.pyx \
    include/Point3D.pyx \
    include/Profile1D.pyx \
    include/Profile2D.pyx \
    include/ProfileBin1D.pyx \
    include/ProfileBin2D.pyx \
    include/Scatter1D.pyx \
    include/Scatter2D.pyx \
    include/Scatter3D.pyx

STATIC_SRC += \
    $(CYTHON_SRC) \
    mktemplates \
    core.cpp core.h util.cpp \
    rootcompat.cpp \
    .made_pyx_templates

BUILT_SRC = \
    include/Bin1D.pxi \
    include/Bin2D.pxi \
    include/Axis1D.pxi \
    include/Axis2D.pxi \
    include/generated/Bin1D_Dbn1D.pyx \
    include/generated/Bin1D_Dbn2D.pyx \
    include/generated/Bin2D_Dbn2D.pyx \
    include/generated/Bin2D_Dbn3D.pyx

MAINTAINERCLEANFILES = \
  core.h core.cpp util.cpp rootcompat.cpp .made_pyx_templates


EXTRA_DIST = $(STATIC_SRC) $(BUILT_SRC)

if WITH_CYTHON

core.cpp core.h: $(CYTHON_SRC) .made_pyx_templates
	$(CYTHON) $(srcdir)/core.pyx --cplus -2 \
       -I $(srcdir) -I $(srcdir)/include \
       -I $(builddir) -I $(builddir)/include \
       -o $@

util.cpp: util.pxd util.pyx .made_pyx_templates
	$(CYTHON) $(srcdir)/util.pyx --cplus -2 \
       -I $(srcdir) -I $(srcdir)/include \
       -I $(builddir) -I $(builddir)/include \
       -o $@

rootcompat.cpp: util.pxd rootcompat.pxd rootcompat.pyx .made_pyx_templates
	$(CYTHON) $(srcdir)/rootcompat.pyx --cplus -2 \
       -I $(srcdir) -I $(srcdir)/include \
       -I $(builddir) -I $(builddir)/include \
       -o $@

else

core.cpp core.h: $(CYTHON_SRC)
	@echo "Not (re)generating core.cpp since Cython is not installed"

util.cpp: util.pxd util.pyx
	@echo "Not (re)generating util.cpp since Cython is not installed"

rootcompat.cpp: util.pxd rootcompat.pxd rootcompat.pyx
	@echo "Not (re)generating rootcompat.cpp since Cython is not installed"

endif

BUILT_SOURCES = util.cpp core.cpp core.h rootcompat.cpp

## Make text-substitution .pyx files from templates -- should only run once
.made_pyx_templates: include/Bin1D_DBN.pyx include/Bin2D_DBN.pyx include/Axis1D_BIN1D_DBN.pyx include/Axis2D_BIN2D_DBN.pyx
	$(srcdir)/mktemplates $^
	touch $@
