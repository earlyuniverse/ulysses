from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map
from libcpp cimport bool
from libcpp.string cimport string
from cython.operator cimport dereference as deref

cdef extern from "YODA/Config/YodaConfig.h" namespace "YODA":
     string version()


# Import the error handling C++ routine
cdef extern from "errors.hh":
    # Have a look in errors.cpp for implementation specifics
    void yodaerr "translate_yoda_error" ()

ctypedef map[string, string] Annotations
ctypedef double (*dbl_dbl_fptr) (double)
ctypedef map[string, pair[double,double]] errMap


# Math utils {{{
cdef extern from "YODA/Utils/MathUtils.h" namespace "YODA":
    # bool isZero(double a, double tolerance)
    # bool fuzzyEquals(double a, double b, double tolerance)
    # bool fuzzyGtrEquals(double a, double b, double tolerance)
    # bool fuzzyLessEquals(double a, double b, double tolerance)
    vector[double] linspace(size_t nbins, double start, double end)
    vector[double] logspace(size_t nbins, double start, double end)
    int index_between(double&, vector[double]& binedges)
    double mean(vector[int]& sample)
    double covariance(vector[int]& sample1, vector[int]& sample2)
    double correlation(vector[int]& sample1, vector[int]& sample2)
# }}}


# Dbn0D {{{
cdef extern from "YODA/Dbn0D.h" namespace "YODA":
    cdef cppclass Dbn0D:
        Dbn0D ()
        Dbn0D (Dbn0D)

        void fill(double weight, double fraction)
        void reset()
        void scaleW(double)

        # Raw distribution running sums
        unsigned long numEntries() except +yodaerr
        double effNumEntries() except +yodaerr
        double sumW() except +yodaerr
        double sumW2() except +yodaerr

        double errW() except +yodaerr
        double relErrW() except +yodaerr

        Dbn0D operator+ (Dbn0D)
        Dbn0D operator- (Dbn0D)
        # TODO: += and -= operators

#}}} Dbn0D


# Dbn1D {{{
cdef extern from "YODA/Dbn1D.h" namespace "YODA":
    cdef cppclass Dbn1D:
        Dbn1D ()
        Dbn1D (Dbn1D)

        void fill(double val, double weight, double fraction)
        void reset()
        void scaleW(double)
        void scaleX(double)

        double errW() except +yodaerr
        double relErrW() except +yodaerr

        double xMean() except +yodaerr
        double xVariance() except +yodaerr
        double xStdDev() except +yodaerr
        double xStdErr() except +yodaerr
        double xRMS() except +yodaerr

        # Raw distribution running sums
        unsigned long numEntries() except +yodaerr
        double effNumEntries() except +yodaerr
        double sumW() except +yodaerr
        double sumW2() except +yodaerr
        double sumWX() except +yodaerr
        double sumWX2() except +yodaerr

        Dbn1D operator+ (Dbn1D)
        Dbn1D operator- (Dbn1D)
        # TODO: += and -= operators

#}}} Dbn1D


# Dbn2D {{{
cdef extern from "YODA/Dbn2D.h" namespace "YODA":
    cdef cppclass Dbn2D:
        Dbn2D ()
        Dbn2D (Dbn2D)

        void fill(double x, double y, double weight, double fraction) except +yodaerr
        void reset() except +yodaerr
        void scaleW(double) except +yodaerr
        void scaleX(double) except +yodaerr
        void scaleY(double) except +yodaerr
        void scaleXY(double, double) except +yodaerr

        double errW() except +yodaerr
        double relErrW() except +yodaerr

        double xMean() except +yodaerr
        double xVariance() except +yodaerr
        double xStdDev() except +yodaerr
        double xStdErr() except +yodaerr
        double xRMS() except +yodaerr

        double yMean() except +yodaerr
        double yVariance() except +yodaerr
        double yStdDev() except +yodaerr
        double yStdErr() except +yodaerr
        double yRMS() except +yodaerr

        # Raw distribution running sums
        unsigned long numEntries() except +yodaerr
        double effNumEntries() except +yodaerr
        double sumW() except +yodaerr
        double sumW2() except +yodaerr
        double sumWX() except +yodaerr
        double sumWX2() except +yodaerr
        double sumWY() except +yodaerr
        double sumWY2() except +yodaerr
        double sumWXY() except +yodaerr

        # Operators
        void flipXY() except +yodaerr
        Dbn1D transformX() except +yodaerr
        Dbn1D transformY() except +yodaerr

        Dbn2D operator + (Dbn2D)
        Dbn2D operator - (Dbn2D)
        # TODO: += and -= operators

#}}} Dbn2D


# Dbn3D {{{
cdef extern from "YODA/Dbn3D.h" namespace "YODA":
    cdef cppclass Dbn3D:
        Dbn3D ()
        Dbn3D (Dbn3D)
        void fill(double x, double y, double z, double weight, double fraction)
        void reset()

        void scaleW(double)
        void scaleX(double)
        void scaleY(double)
        void scaleZ(double)
        # void scaleXY(double, double)
        # void scaleYZ(double, double)
        # void scaleXZ(double, double)
        void scaleXYZ(double, double, double)

        double errW() except +yodaerr
        double relErrW() except +yodaerr

        double xMean()
        double xVariance()
        double xStdDev()
        double xStdErr()
        double xRMS()

        double yMean()
        double yVariance()
        double yStdDev()
        double yStdErr()
        double yRMS()

        double zMean()
        double zVariance()
        double zStdDev()
        double zStdErr()
        double zRMS()

        # Raw distribution running sums
        unsigned long numEntries()
        double effNumEntries()
        double sumW()
        double sumW2()

        double sumWX()
        double sumWX2()

        double sumWY()
        double sumWY2()

        double sumWZ()
        double sumWZ2()

        double sumWXY()
        double sumWXZ()
        double sumWYZ()

        double sumWXYZ()

        # Operators
        void flipXY()
        void flipXZ()
        void flipYZ()

        Dbn1D transformX()
        Dbn1D transformY()
        Dbn1D transformZ()

        Dbn3D operator + (Dbn3D)
        Dbn3D operator - (Dbn3D)
        # TODO: += and -= operators

#}}} Dbn3D


# Point {{{

cdef extern from "YODA/Point.h" namespace "YODA":
    cdef cppclass Point:

        int dim() except +yodaerr

        double val(size_t i) except +yodaerr
        void setVal(size_t i, double val) except +yodaerr

        pair[double,double] errs(size_t i) except +yodaerr
        pair[double,double] errs(size_t i, string source) except +yodaerr
        double errMinus(size_t i) except +yodaerr
        double errMinus(size_t i, string source) except +yodaerr
        void setErrMinus(size_t i, double eminus) except +yodaerr
        void setErrMinus(size_t i, double eminus, string source) except +yodaerr
        double errPlus(size_t i) except +yodaerr
        double errPlus(size_t i, string source) except +yodaerr
        void setErrPlus(size_t i, double eplus) except +yodaerr
        void setErrPlus(size_t i, double eplus, string source) except +yodaerr
        double errAvg(size_t i) except +yodaerr
        double errAvg(size_t i, string source) except +yodaerr

        void setErr(size_t i, double e) except +yodaerr
        void setErr(size_t i, double e, string source) except +yodaerr
        # void setErrs(size_t i, double e) except +yodaerr
        # void setErrs(size_t i, double eminus, double eplus) except +yodaerr
        void setErrs(size_t i, pair[double,double]& e) except +yodaerr
        void setErrs(size_t i, pair[double,double]& e, string source) except +yodaerr

        # void set(size_t i, double val, double e) except +yodaerr
        # void set(size_t i, double val, double eminus, double eplus) except +yodaerr
        void set(size_t i, double val, pair[double,double]& e) except +yodaerr
        void set(size_t i, double val, pair[double,double]& e, string source) except +yodaerr
        void scale(size_t i, double scale) except +yodaerr

        errMap errMap() except +yodaerr
        void rmVariations() except +yodaerr

#}}} Point


# Point1D {{{
cdef extern from "YODA/Point1D.h" namespace "YODA":
    cdef cppclass Point1D(Point):
        Point1D () except +yodaerr
        Point1D (Point1D p) except +yodaerr
        Point1D (double x, double exminus, double explus) except +yodaerr
        Point1D (double x, double exminus, double explus, string source) except +yodaerr

        double x() except +yodaerr
        void setX(double x) except +yodaerr

        pair[double,double] xErrs() except +yodaerr
        pair[double,double] xErrs(string source) except +yodaerr
        void setXErrs(pair[double, double]&) except +yodaerr
        void setXErrs(pair[double, double]&, string source) except +yodaerr
        double xErrAvg() except +yodaerr
        double xErrAvg(string source) except +yodaerr

        double xMin() except +yodaerr
        double xMin(string source) except +yodaerr
        double xMax() except +yodaerr
        double xMax(string source) except +yodaerr

        void scaleX(double) except +yodaerr
        void scale(size_t i, double scale) except +yodaerr
        void updateTotalUncertainty() except +yodaerr

        bool operator == (Point1D) except +yodaerr
        bool operator != (Point1D b) except +yodaerr
        bool operator < (Point1D b) except +yodaerr
        bool operator <= (Point1D b) except +yodaerr
        bool operator > (Point1D b) except +yodaerr
        bool operator >= (Point1D b) except +yodaerr
# }}} Point1D


# Point2D {{{
cdef extern from "YODA/Point2D.h" namespace "YODA":
    cdef cppclass Point2D(Point):
        Point2D () except +yodaerr
        Point2D (Point2D p) except +yodaerr
        Point2D (double x, double y,
                 double exminus, double explus,
                 double eyminus, double eyplus) except +yodaerr
        Point2D (double x, double y,
                 double exminus, double explus,
                 double eyminus, double eyplus, string source) except +yodaerr

        double x() except +yodaerr
        double y() except +yodaerr
        void setX(double x) except +yodaerr
        void setY(double y) except +yodaerr
        pair[double,double] xy() except +yodaerr
        void setXY(pair[double,double]&) except +yodaerr

        pair[double,double] xErrs() except +yodaerr
        pair[double,double] yErrs() except +yodaerr
        pair[double,double] yErrs(string source) except +yodaerr
        void setXErrs(pair[double, double]&) except +yodaerr
        void setYErrs(pair[double, double]&) except +yodaerr
        void setYErrs(pair[double, double]&, string source) except +yodaerr
        double xErrAvg() except +yodaerr
        double yErrAvg() except +yodaerr
        double yErrAvg(string source) except +yodaerr

        double xMin() except +yodaerr
        double xMax() except +yodaerr
        double yMin() except +yodaerr
        double yMin(string source) except +yodaerr
        double yMax() except +yodaerr
        double yMax(string source) except +yodaerr

        void scaleX(double) except +yodaerr
        void scaleY(double) except +yodaerr
        void scaleXY(double, double) except +yodaerr
        void scale(size_t i, double scale) except +yodaerr

        void updateTotalUncertainty() except +yodaerr

        bool operator == (Point2D) except +yodaerr
        bool operator != (Point2D b) except +yodaerr
        bool operator < (Point2D b) except +yodaerr
        bool operator <= (Point2D b) except +yodaerr
        bool operator > (Point2D b) except +yodaerr
        bool operator >= (Point2D b) except +yodaerr
# }}} Point2D


# Point3D {{{
cdef extern from "YODA/Point3D.h" namespace "YODA":
    cdef cppclass Point3D(Point):
        Point3D () except +yodaerr
        Point3D (Point3D& p) except +yodaerr
        Point3D (double x, double y, double z,
                 double exminus, double explus,
                 double eyminus, double eyplus,
                 double ezminus, double ezplus) except +yodaerr
        Point3D (double x, double y, double z,
                 double exminus, double explus,
                 double eyminus, double eyplus,
                 double ezminus, double ezplus, string source) except +yodaerr

        double x() except +yodaerr
        double y() except +yodaerr
        double z() except +yodaerr
        void setX(double x) except +yodaerr
        void setY(double y) except +yodaerr
        void setZ(double z) except +yodaerr

        pair[double,double] xErrs() except +yodaerr
        pair[double,double] yErrs() except +yodaerr
        pair[double,double] zErrs() except +yodaerr
        pair[double,double] zErrs(string source) except +yodaerr
        void setXErrs(pair[double, double]&) except +yodaerr
        void setYErrs(pair[double, double]&) except +yodaerr
        void setZErrs(pair[double, double]&) except +yodaerr
        void setZErrs(pair[double, double]&, string source) except +yodaerr
        double xErrAvg()
        double yErrAvg()
        double zErrAvg()
        double zErrAvg(string source)

        double xMin() except +yodaerr
        double xMax() except +yodaerr
        double yMin() except +yodaerr
        double yMax() except +yodaerr
        double zMin() except +yodaerr
        double zMin(string source) except +yodaerr
        double zMax() except +yodaerr
        double zMax(string source) except +yodaerr

        void scaleX(double) except +yodaerr
        void scaleY(double) except +yodaerr
        void scaleZ(double) except +yodaerr
        void scaleXYZ(double, double, double) except +yodaerr
        void scale(size_t i, double scale) except +yodaerr

        void updateTotalUncertainty() except +yodaerr

        bool operator == (Point3D b)
        bool operator != (Point3D b)
        bool operator < (Point3D b)
        bool operator <= (Point3D b)
        bool operator > (Point3D b)
        bool operator >= (Point3D b)
#}}} Point3D





# Bin {{{
cdef extern from "YODA/Bin.h" namespace "YODA":
    cdef cppclass Bin:
        int dim() except +yodaerr
        int fillDim() except +yodaerr
        unsigned long numEntries() except +yodaerr
        double effNumEntries() except +yodaerr
        double sumW() except +yodaerr
        double sumW2() except +yodaerr
# }}} Bin



#Bin1D {{{
cdef extern from "YODA/Bin1D.h" namespace "YODA":
    cdef cppclass Bin1D[DBN](Bin):
        Bin1D(pair[double, double] edges) except +yodaerr
        Bin1D(pair[double, double] edges, DBN dbn) except +yodaerr
        Bin1D(Bin1D) except +yodaerr

        # THIS IS A CYTHON LIMITATION... DO NOT CALL THIS
        Bin1D() # (DO NOT CALL THIS DO NOT CALL THIS) ###
        #################################################

        #We're fine as long as we don't try to instantiate these from Python

        # void scaleW(double scale) except +yodaerr
        # void scaleX(double scale) except +yodaerr
        void reset()  except +yodaerr

        pair[double, double] edges() except +yodaerr

        double xMin() except +yodaerr
        double xMax() except +yodaerr
        double xMid() except +yodaerr
        double xWidth() except +yodaerr

        double xFocus() except +yodaerr

        # x statistics
        double xMean() except +yodaerr
        double xVariance() except +yodaerr
        double xStdDev() except +yodaerr
        double xStdErr() except +yodaerr
        double xRMS() except +yodaerr

        # raw statistics
        double sumWX() except +yodaerr
        double sumWX2() except +yodaerr

        void merge (Bin1D&) except +yodaerr
        Bin1D operator + (Bin1D&)
        Bin1D operator - (Bin1D&)

ctypedef Bin1D[Dbn1D] Bin1D_Dbn1D
ctypedef Bin1D[Dbn2D] Bin1D_Dbn2D
ctypedef Bin1D[Dbn3D] Bin1D_Dbn3D
#}}} Bin1D


# Bin2D {{{
cdef extern from "YODA/Bin2D.h" namespace "YODA":
    cdef cppclass Bin2D[DBN](Bin):
        Bin2D(pair[double, double] xedges, pair[double, double] yedges) except+
        Bin2D(Bin2D bin) except +yodaerr

        # CYTHON HACK DO NOT CALL THIS IT DOES NOT EXIST
        Bin2D() # (DO NOT CALL DO NOT CALL)
        ################################################

        # void scaleW(double scale) except +yodaerr
        # void scaleXY(double, double) except +yodaerr
        void reset()  except +yodaerr

        pair[double, double] xEdges() except +yodaerr
        pair[double, double] yEdges() except +yodaerr

        double xMin() except +yodaerr
        double yMin() except +yodaerr
        double xMax() except +yodaerr
        double yMax() except +yodaerr
        double xMid() except +yodaerr
        double yMid() except +yodaerr
        double xWidth() except +yodaerr
        double yWidth() except +yodaerr
        double area() except +yodaerr

        double xFocus() except +yodaerr
        double yFocus() except +yodaerr

        pair[double, double] xyFocus() except +yodaerr
        pair[double, double] xyMid() except +yodaerr

        # x statistics
        double xMean() except +yodaerr
        double xVariance() except +yodaerr
        double xStdDev() except +yodaerr
        double xStdErr() except +yodaerr
        double xRMS() except +yodaerr

        double yMean() except +yodaerr
        double yVariance() except +yodaerr
        double yStdDev() except +yodaerr
        double yStdErr() except +yodaerr
        double yRMS() except +yodaerr

        # Raw statistics
        double sumWX() except +yodaerr
        double sumWY() except +yodaerr
        double sumWXY() except +yodaerr
        double sumWX2() except +yodaerr
        double sumWY2() except +yodaerr

        #void merge(Bin2D) except +yodaerr
        Bin2D operator + (Bin2D)
        Bin2D operator - (Bin2D)

        int adjacentTo(Bin2D) except +yodaerr

ctypedef Bin2D[Dbn2D] Bin2D_Dbn2D
ctypedef Bin2D[Dbn3D] Bin2D_Dbn3D
# }}} Bin2D



# HistoBin1D {{{
cdef extern from "YODA/HistoBin1D.h" namespace "YODA":
    cdef cppclass HistoBin1D(Bin1D_Dbn1D):
        HistoBin1D(double lowedge, double highedge) except +yodaerr
        HistoBin1D(HistoBin1D) except +yodaerr
        # void fill(double x, double weight, double fraction) except +yodaerr
        # void fillBin(double weight, double fraction) except +yodaerr

        double area() except +yodaerr
        double height() except +yodaerr
        double areaErr() except +yodaerr
        double heightErr() except +yodaerr
        double relErr() except +yodaerr

        HistoBin1D operator+(HistoBin1D)
        HistoBin1D operator-(HistoBin1D)

#}}} HistoBin1D

cdef extern from "merge.hh":
    void HistoBin1D_iadd_HistoBin1D "cython_iadd" (HistoBin1D*, HistoBin1D*)
    void HistoBin1D_isub_HistoBin1D "cython_isub" (HistoBin1D*, HistoBin1D*)
    # void HistoBin1D_imul_dbl "cython_imul_dbl" (HistoBin1D*, double)
    # void HistoBin1D_idiv_dbl "cython_idiv_dbl" (HistoBin1D*, double)
    HistoBin1D* HistoBin1D_add_HistoBin1D "cython_add" (HistoBin1D*, HistoBin1D*)
    HistoBin1D* HistoBin1D_sub_HistoBin1D "cython_sub" (HistoBin1D*, HistoBin1D*)
    HistoBin1D* HistoBin1D_div_HistoBin1D "cython_div" (HistoBin1D*, HistoBin1D*)


# HistoBin2D {{{
cdef extern from "YODA/HistoBin2D.h" namespace "YODA":
    cdef cppclass HistoBin2D(Bin2D_Dbn2D):
        HistoBin2D(double xmin, double xmax, double ymin, double ymax) except +yodaerr
        HistoBin2D(HistoBin2D) except +yodaerr

        # void fill(double x, double y, double weight, double fraction) except +yodaerr
        # void fillBin(double weight, double fraction) except +yodaerr
        void reset()

        # Accessors
        double volume() except +yodaerr
        double volumeErr() except +yodaerr
        double height() except +yodaerr
        double heightErr() except +yodaerr
        double relErr() except +yodaerr

        HistoBin2D operator+(HistoBin2D)
        HistoBin2D operator-(HistoBin2D)

        #Bin2D_Dbn2D merge(HistoBin2D b)
#}}} HistoBin2D



# ProfileBin1D {{{
cdef extern from "YODA/ProfileBin1D.h" namespace "YODA":

    cdef cppclass ProfileBin1D(Bin1D_Dbn2D):
        ProfileBin1D(ProfileBin1D) except +yodaerr
        ProfileBin1D(double, double) except +yodaerr
        #void fill(double x, double y, double weight, double fraction) except +yodaerr
        #void fillBin(double y, double weight, double fraction) except +yodaerr
        void reset() except +yodaerr

        double mean() except +yodaerr
        double stdDev() except +yodaerr
        double variance() except +yodaerr
        double stdErr() except +yodaerr
        double rms() except +yodaerr

        double sumWY() except +yodaerr
        double sumWY2() except +yodaerr
        ProfileBin1D operator + (ProfileBin1D)
        ProfileBin1D operator - (ProfileBin1D)

        # void scaleY(double) except +yodaerr

# }}} ProfileBin1D

cdef extern from "merge.hh":
    void ProfileBin1D_iadd_ProfileBin1D "cython_iadd" (ProfileBin1D*, ProfileBin1D*)
    void ProfileBin1D_isub_ProfileBin1D "cython_isub" (ProfileBin1D*, ProfileBin1D*)
    # void ProfileBin1D_imul_dbl "cython_imul_dbl" (ProfileBin1D*, double)
    # void ProfileBin1D_idiv_dbl "cython_idiv_dbl" (ProfileBin1D*, double)
    ProfileBin1D* ProfileBin1D_add_ProfileBin1D "cython_add" (ProfileBin1D*, ProfileBin1D*)
    ProfileBin1D* ProfileBin1D_sub_ProfileBin1D "cython_sub" (ProfileBin1D*, ProfileBin1D*)
    ProfileBin1D* ProfileBin1D_div_ProfileBin1D "cython_div" (ProfileBin1D*, ProfileBin1D*)


# ProfileBin2D {{{
cdef extern from "YODA/ProfileBin2D.h" namespace "YODA":

    cdef cppclass ProfileBin2D(Bin2D_Dbn3D):
        ProfileBin2D (ProfileBin2D h) except +yodaerr
        ProfileBin2D (double, double, double, double) except +yodaerr
        # void fill(double x, double y, double z, double weight, double fraction) except +yodaerr
        # void fillBin(double z, double weight, double fraction) except +yodaerr

        double mean() except +yodaerr
        double stdDev() except +yodaerr
        double variance() except +yodaerr
        double stdErr() except +yodaerr
        double rms() except +yodaerr

        double sumWZ() except +yodaerr
        double sumWZ2() except +yodaerr
        ProfileBin2D operator + (ProfileBin2D)
        ProfileBin2D operator - (ProfileBin2D)

        # void scaleZ(double) except +yodaerr

# }}} ProfileBin2D




# AnalysisObject {{{
cdef extern from "YODA/AnalysisObject.h" namespace "YODA":
    cdef cppclass AnalysisObject:
        # Constructors
        AnalysisObject(string type, string path, string title) except +yodaerr
        AnalysisObject(string type, string path, AnalysisObject ao, string title) except +yodaerr
        AnalysisObject()
        #AnalysisObject* newclone() except +yodaerr

        ## String used in automatic type determination
        string type() except +yodaerr

        ## Data object fill- or plot-space dimension
        int dim() except +yodaerr

        ## Annotations
        vector[string] annotations() except +yodaerr
        bool hasAnnotation(string key) except +yodaerr
        string annotation(string key) except +yodaerr
        string annotation(string key, string default) except +yodaerr
        void setAnnotation(string, string) except +yodaerr
        void rmAnnotation(string name) except +yodaerr
        void clearAnnotations() except +yodaerr

        ## Standard annotations
        string title() except +yodaerr
        void setTitle(string title) except +yodaerr
        string path() except +yodaerr
        void setPath(string title) except +yodaerr
        string name() except +yodaerr
# }}} AnalysisObject


# Fillable#{{{
cdef extern from "YODA/Fillable.h" namespace "YODA":
    cdef cppclass Fillable:
        pass
#}}}

# Binned#{{{
cdef extern from "YODA/Binned.h" namespace "YODA":
    cdef cppclass Binned:
        pass
#}}}

# Scatter#{{{
cdef extern from "YODA/Scatter.h" namespace "YODA":
    cdef cppclass Scatter:
        pass
#}}}


cdef extern from "YODA/Utils/sortedvector.h" namespace "YODA::Utils":
    cdef cppclass sortedvector[T](vector):
        sortedvector(vector[T]) except +yodaerr
        void insert(T) except +yodaerr

# TODO: forward declarations for bin-copying constructors


# Counter {{{
cdef extern from "YODA/Counter.h" namespace "YODA":
    cdef cppclass Counter(AnalysisObject,Fillable):
        Counter() except +yodaerr

        Counter(string path, string title) except +yodaerr

        #Counter(Dbn0D dbn, string path, string title) except +yodaerr

        Counter(Counter c, string path)

        Counter clone() except +yodaerr
        Counter* newclone() except +yodaerr


        void reset() except +yodaerr

        void fill(double weight, double fraction) except +yodaerr

        unsigned long numEntries() except +yodaerr
        double effNumEntries() except +yodaerr

        double sumW() except +yodaerr
        double sumW2() except +yodaerr

        double val() except +yodaerr
        double err() except +yodaerr
        double relErr() except +yodaerr

        void scaleW(double) except +yodaerr

        # operator += (Counter)
        # operator -= (Counter)

    Scatter1D Counter_div_Counter "divide" (const Counter&, const Counter&) except +yodaerr
    Scatter1D Counter_eff_Counter "efficiency" (const Counter&, const Counter&) except +yodaerr


cdef extern from "merge.hh":
    void Counter_iadd_Counter "cython_iadd" (Counter*, Counter*)
    void Counter_isub_Counter "cython_isub" (Counter*, Counter*)
    # void Counter_imul_dbl "cython_imul_dbl" (Counter*, double)
    # void Counter_idiv_dbl "cython_idiv_dbl" (Counter*, double)
    Counter* Counter_add_Counter "cython_add" (Counter*, Counter*)
    Counter* Counter_sub_Counter "cython_sub" (Counter*, Counter*)
    #Counter* Counter_div_Counter "cython_div" (Counter*, Counter*)

cdef extern from "YODA/Scatter1D.h" namespace "YODA":
    Scatter1D mkScatter_Counter "YODA::mkScatter" (const Counter&) except +yodaerr

#}}} Counter


# Scatter1D {{{
cdef extern from "YODA/Scatter1D.h" namespace "YODA":
    cdef cppclass Scatter1D(AnalysisObject,Scatter):
        Scatter1D() except +yodaerr

        Scatter1D(string path, string title) except +yodaerr

        Scatter1D(sortedvector[Point1D],
                  string path,
                  string title) except +yodaerr

        Scatter1D(vector[double], vector[double],
                  vector[pair[double, double]],
                  vector[pair[double, double]]) except +yodaerr

        Scatter1D(Scatter1D p, string path)

        Scatter1D clone() except +yodaerr
        Scatter1D* newclone() except +yodaerr


        void reset() except +yodaerr


        size_t numPoints() except +yodaerr
        # TODO: have to ignore exception handling on ref-returning methods until Cython bug is fixed
        vector[Point1D]& points() #except +yodaerr
        Point1D& point(size_t index) #except +yodaerr

        void addPoint(Point1D&) #except +yodaerr
        void addPoint(double) #except +yodaerr
        void addPoint(double, const pair[double, double]&) #except +yodaerr

        void addPoints(const sortedvector[Point1D]&) #except +yodaerr

        void rmPoint(size_t) #except +yodaerr
        void rmPoints(const vector[size_t]) #except +yodaerr

        void combineWith(const Scatter1D&) #except +yodaerr
        void combineWith(const vector[Scatter1D]&) #except +yodaerr

        void scaleX(double) except +yodaerr
        void scale(size_t i, double scale) except +yodaerr

        void parseVariations() except +yodaerr
        void writeVariationsToAnnotations() except +yodaerr
        void updateTotalUncertainty() except +yodaerr
        void rmVariations() except +yodaerr
        vector[string] variations() except +yodaerr

        vector[vector[double]] covarianceMatrix(bool) except +yodaerr

    void Scatter1D_transformX "YODA::transformX" (Scatter1D&, dbl_dbl_fptr)

#}}} Scatter1D

# cdef extern from "merge.hh":
#     Scatter2D* Scatter2D_add_Scatter2D "cython_add" (Scatter2D*, Scatter2D*)
#     Scatter2D* Scatter2D_sub_Scatter2D "cython_sub" (Scatter2D*, Scatter2D*)

cdef extern from "YODA/Scatter1D.h" namespace "YODA":
    Scatter1D mkScatter_Scatter1D "YODA::mkScatter" (const Scatter1D&) except +yodaerr


# Scatter2D {{{
cdef extern from "YODA/Scatter2D.h" namespace "YODA":
    cdef cppclass Scatter2D(AnalysisObject,Scatter):
        Scatter2D() except +yodaerr

        Scatter2D(string path, string title) except +yodaerr

        Scatter2D(sortedvector[Point2D],
                  string path,
                  string title) except +yodaerr

        Scatter2D(vector[double], vector[double],
                  vector[pair[double, double]],
                  vector[pair[double, double]]) except +yodaerr

        Scatter2D(Scatter2D p, string path)

        Scatter2D clone() except +yodaerr
        Scatter2D* newclone() except +yodaerr


        void reset() except +yodaerr


        size_t numPoints() except +yodaerr
        # TODO: have to ignore exception handling on ref-returning methods until Cython bug is fixed
        vector[Point2D]& points() #except +yodaerr
        Point2D& point(size_t index) #except +yodaerr

        void addPoint(Point2D&) #except +yodaerr
        void addPoint(double, double) #except +yodaerr
        void addPoint(double, double,
                      const pair[double, double]&, const pair[double, double]&) #except +yodaerr

        void addPoints(const sortedvector[Point2D]&) #except +yodaerr

        void rmPoint(size_t) #except +yodaerr
        void rmPoints(const vector[size_t]) #except +yodaerr

        void combineWith(const Scatter2D&) #except +yodaerr
        void combineWith(const vector[Scatter2D]&) #except +yodaerr

        void scaleX(double) except +yodaerr
        void scaleY(double) except +yodaerr
        void scaleXY(double, double) except +yodaerr
        void scale(size_t i, double scale) except +yodaerr

        void writeVariationsToAnnotations() except +yodaerr
        void updateTotalUncertainty() except +yodaerr
        void parseVariations() except +yodaerr
        void rmVariations() except +yodaerr
        vector[string] variations() except +yodaerr

        vector[vector[double]] covarianceMatrix(bool) except +yodaerr

    void Scatter2D_transformX "YODA::transformX" (Scatter2D&, dbl_dbl_fptr)
    void Scatter2D_transformY "YODA::transformY" (Scatter2D&, dbl_dbl_fptr)

#}}} Scatter2D

# cdef extern from "merge.hh":
#     Scatter2D* Scatter2D_add_Scatter2D "cython_add" (Scatter2D*, Scatter2D*)
#     Scatter2D* Scatter2D_sub_Scatter2D "cython_sub" (Scatter2D*, Scatter2D*)

cdef extern from "YODA/Scatter2D.h" namespace "YODA":
    Scatter2D mkScatter_Scatter2D "YODA::mkScatter" (const Scatter2D&) except +yodaerr



# Scatter3D {{{
cdef extern from "YODA/Scatter3D.h" namespace "YODA":
    cdef cppclass Scatter3D(AnalysisObject,Scatter):
        Scatter3D() except +yodaerr

        Scatter3D(string path, string title) except +yodaerr

        Scatter3D(sortedvector[Point3D],
                  string path,
                  string title) except +yodaerr

        Scatter3D(vector[double], vector[double],
                  vector[pair[double, double]],
                  vector[pair[double, double]],
                  vector[pair[double, double]]) except +yodaerr

        Scatter3D(Scatter3D p, string path)

        Scatter3D clone() except +yodaerr
        Scatter3D* newclone() except +yodaerr


        void reset() except +yodaerr


        size_t numPoints() except +yodaerr
        # TODO: have to ignore exception handling on ref-returning methods until Cython bug is fixed
        sortedvector[Point3D]& points() #except +yodaerr
        Point3D& point(size_t index) #except +yodaerr

        void addPoint(Point3D&) #except +yodaerr
        void addPoint(double, double, double) #except +yodaerr
        void addPoint(double, double, double,
                      const pair[double, double]&, const pair[double, double]&, const pair[double, double]&) #except +yodaerr

        void addPoints(const sortedvector[Point3D]&) #except +yodaerr

        void rmPoint(size_t) #except +yodaerr
        void rmPoints(const vector[size_t]) #except +yodaerr

        void combineWith(const Scatter3D&) #except +yodaerr
        void combineWith(const vector[Scatter3D]&) #except +yodaerr

        void scaleX(double) except +yodaerr
        void scaleY(double) except +yodaerr
        void scaleZ(double) except +yodaerr
        void scaleXYZ(double, double, double) except +yodaerr
        void scale(size_t i, double scale) except +yodaerr

        void parseVariations() except +yodaerr
        void updateTotalUncertainty() except +yodaerr
        void writeVariationsToAnnotations() except +yodaerr
        void rmVariations() except +yodaerr
        vector[string] variations() except +yodaerr

        vector[vector[double]] covarianceMatrix() except +yodaerr
        vector[vector[double]] covarianceMatrix(bool) except +yodaerr

    void Scatter3D_transformX "YODA::transformX" (Scatter3D&, dbl_dbl_fptr)
    void Scatter3D_transformY "YODA::transformY" (Scatter3D&, dbl_dbl_fptr)
    void Scatter3D_transformZ "YODA::transformZ" (Scatter3D&, dbl_dbl_fptr)

#}}} Scatter3D

# cdef extern from "merge.hh":
#     Scatter3D* Scatter3D_add_Scatter3D "cython_add" (Scatter3D*, Scatter3D*)
#     Scatter3D* Scatter3D_sub_Scatter3D "cython_sub" (Scatter3D*, Scatter3D*)

cdef extern from "YODA/Scatter3D.h" namespace "YODA":
    Scatter3D mkScatter_Scatter3D "YODA::mkScatter" (const Scatter3D&) except +yodaerr


# Histo1D#{{{
cdef extern from "YODA/Histo1D.h" namespace "YODA":
    cdef cppclass Histo1D(AnalysisObject,Fillable,Binned):
        Histo1D() except +yodaerr

        Histo1D(string path, string title) except +yodaerr

        Histo1D(size_t nbins,
                double lower,
                double upper,
                string path,
                string title) except +yodaerr

        Histo1D(vector[double] binedges,
                string path,
                string title) except +yodaerr

        Histo1D(vector[Bin] bins, string path, string title) except +yodaerr

        Histo1D(Histo1D h, string path) except +yodaerr

        #Histo1D(Profile1D p, string path)

        #Histo1D(Scatter2D p, string path)

        Histo1D clone() except +yodaerr
        Histo1D* newclone() except +yodaerr


        void reset() except +yodaerr

        void fill(double x, double weight, double fraction) except +yodaerr
        void fillBin(size_t i, double weight, double fraction) except +yodaerr

        void scaleW(double s) except +yodaerr
        void normalize(double normto, bool includeoverflows) except +yodaerr

        void mergeBins(size_t, size_t) except +yodaerr
        void rebinBy(unsigned int n, size_t begin, size_t end) except +yodaerr
        void rebinTo(vector[double] edges) except +yodaerr

        void addBin(double, double) except +yodaerr
        void addBins(vector[double] edges) except +yodaerr
        void rmBin(size_t index) except +yodaerr

        vector[double] xEdges() except +yodaerr
        vector[double] xWidths() except +yodaerr

        double xMin() except +yodaerr
        double xMax() except +yodaerr

        size_t numBins() except +yodaerr
        size_t numBinsX() except +yodaerr

        vector[HistoBin1D]& bins()
        int binIndexAt(double x) except +yodaerr
        const HistoBin1D& bin(size_t ix)
        const HistoBin1D& binAt(double x) except +yodaerr

        # TODO: Some Cython mapping problem?
        Dbn1D& totalDbn()
        Dbn1D& underflow()
        Dbn1D& overflow()

        # Whole histo data
        double integral(bool)
        double integralTo(int, bool)
        double integralRange(int, int)

        unsigned long numEntries(bool)
        double effNumEntries(bool)
        double sumW(bool)
        double sumW2(bool)

        double xMean(bool)
        double xVariance(bool)
        double xStdDev(bool)
        double xStdErr(bool)
        double xRMS(bool)

        # operator == (Histo1D)
        # operator != (Histo1D)
        operator + (Histo1D)
        operator - (Histo1D)
        operator / (Histo1D)

    Scatter2D Histo1D_toIntegral "toIntegralHisto" (const Histo1D& h, bool includeunderflow) except +yodaerr
    Scatter2D Histo1D_toIntegralEff "toIntegralEfficiencyHisto" (const Histo1D& h, bool includeunderflow, bool includeoverflow) except +yodaerr
    Scatter2D Histo1D_div_Histo1D "divide" (const Histo1D&, const Histo1D&) except +yodaerr
    Scatter2D Histo1D_eff_Histo1D "efficiency" (const Histo1D&, const Histo1D&) except +yodaerr

cdef extern from "merge.hh":
    void Histo1D_iadd_Histo1D "cython_iadd" (Histo1D*, Histo1D*)
    void Histo1D_isub_Histo1D "cython_isub" (Histo1D*, Histo1D*)
    # void Histo1D_imul_dbl "cython_imul_dbl" (Histo1D*, double)
    # void Histo1D_idiv_dbl "cython_idiv_dbl" (Histo1D*, double)
    Histo1D* Histo1D_add_Histo1D "cython_add" (Histo1D*, Histo1D*)
    Histo1D* Histo1D_sub_Histo1D "cython_sub" (Histo1D*, Histo1D*)
    Histo1D* Histo1D_div_Histo1D "cython_div" (Histo1D*, Histo1D*)

cdef extern from "YODA/Scatter2D.h" namespace "YODA":
    Scatter2D mkScatter_Histo1D "YODA::mkScatter" (const Histo1D&, bool, bool, double, double) except +yodaerr

#}}} Histo1D



# Histo2D {{{
cdef extern from "YODA/Histo2D.h" namespace "YODA":
    cdef cppclass Histo2D(AnalysisObject,Fillable,Binned):
        Histo2D() except +yodaerr

        Histo2D(string path, string title) except +yodaerr

        Histo2D(size_t nBinsX, double lowerX, double upperX,
                size_t nBinsY, double lowerY, double upperY,
                string path, string title) except +yodaerr

        Histo2D(vector[double] xedges, vector[double] yedges,
                string path, string title) except +yodaerr

        Histo2D(Histo2D, string path)
        #Histo2D(Profile1D p, string path)
        #Histo2D(Scatter2D p, string path)

        Histo2D clone() except +yodaerr
        Histo2D* newclone() except +yodaerr


        # TODO: add missing functions and enable refs + exceptions when Cython allows

        void reset() except +yodaerr

        void fill(double x, double y, double weight, double fraction) except +yodaerr
        void fillBin(size_t i, double weight, double fraction) except +yodaerr

        void normalize(double normto, bool includeoverflows) except +yodaerr

        void scaleW(double scalefactor) except +yodaerr
        void scaleXY(double, double)

        # void mergeBins(size_t, size_t) except +yodaerr
        # void rebin(unsigned int n) except +yodaerr

        size_t numBins() except +yodaerr
        size_t numBinsX() except +yodaerr
        size_t numBinsY() except +yodaerr

        vector[HistoBin2D]& bins() #except +yodaerr
        int binIndexAt(double x, double y) except +yodaerr
        const HistoBin2D& bin(size_t ix) #except +yodaerr
        const HistoBin2D& binAt(double x, double y) #except +yodaerr

        void addBin(const pair[double, double]&, const pair[double, double]&)
        void addBins(const vector[HistoBin2D]&)
        void addBin(double, double) except +yodaerr
        void addBins(const vector[double]& edges) except +yodaerr
        void rmBin(size_t index) except +yodaerr

        vector[double] xEdges() except +yodaerr
        vector[double] yEdges() except +yodaerr
        vector[double] xWidths() except +yodaerr
        vector[double] yWidths() except +yodaerr

        double xMin() except +yodaerr
        double xMax() except +yodaerr
        double yMin() except +yodaerr
        double yMax() except +yodaerr

        # Dbn2D& outflow(int, int) #except +yodaerr

        # Whole histo data
        Dbn2D& totalDbn() #except +yodaerr
        double integral(bool)
        unsigned long numEntries(bool)
        double effNumEntries(bool)
        double sumW(bool)
        double sumW2(bool)

        double xMean(bool)
        double yMean(bool)
        double xVariance(bool)
        double yVariance(bool)
        double xStdDev(bool)
        double yStdDev(bool)
        double xStdErr(bool)
        double yStdErr(bool)
        double xRMS(bool)
        double yRMS(bool)

        # operator == (Histo2D)
        # operator != (Histo2D)
        operator + (Histo2D)
        operator - (Histo2D)
        operator / (Histo2D)

    Scatter3D Histo2D_div_Histo2D "divide" (const Histo2D&, const Histo2D&) except +yodaerr
    Scatter3D Histo2D_eff_Histo2D "efficiency" (const Histo2D&, const Histo2D&) except +yodaerr

cdef extern from "merge.hh":
    void Histo2D_iadd_Histo2D "cython_iadd" (Histo2D*, Histo2D*)
    void Histo2D_isub_Histo2D "cython_isub" (Histo2D*, Histo2D*)
    # void Histo2D_imul_dbl "cython_imul_dbl" (Histo2D*, double)
    # void Histo2D_idiv_dbl "cython_idiv_dbl" (Histo2D*, double)
    Histo2D* Histo2D_add_Histo2D "cython_add" (Histo2D*, Histo2D*)
    Histo2D* Histo2D_sub_Histo2D "cython_sub" (Histo2D*, Histo2D*)
    Histo2D* Histo2D_div_Histo2D "cython_div" (Histo2D*, Histo2D*)

cdef extern from "YODA/Scatter3D.h" namespace "YODA":
    Scatter3D mkScatter_Histo2D "YODA::mkScatter" (const Histo2D&, bool, bool) except +yodaerr

# Histo2D }}}





# Profile1D {{{
cdef extern from "YODA/Profile1D.h" namespace "YODA":
    cdef cppclass Profile1D(AnalysisObject,Fillable,Binned):
        Profile1D() except +yodaerr

        Profile1D(string path, string title) except +yodaerr

        Profile1D(size_t nxbins,
                double xlower,
                double xupper,
                string path,
                string title) except +yodaerr

        Profile1D(vector[double] xbinedges,
                string path,
                string title) except +yodaerr

        Profile1D(Profile1D p, string path) except +yodaerr

        Profile1D(Scatter2D s, string path) except +yodaerr

        #Profile1D(Histo1D p, string path)

        Profile1D clone() except +yodaerr
        Profile1D* newclone() except +yodaerr


        void reset() except +yodaerr

        void fill(double x, double y, double weight, double fraction) except +yodaerr
        void fillBin(size_t i, double y, double weight, double fraction) except +yodaerr

        void scaleW(double s) except +yodaerr
        void scaleY(double s) except +yodaerr

        void mergeBins(size_t, size_t) except +yodaerr
        void rebinBy(unsigned int n, size_t begin, size_t end) except +yodaerr
        void rebinTo(vector[double] edges) except +yodaerr

        void addBin(double, double) except +yodaerr
        void addBins(vector[double] edges) except +yodaerr
        void rmBin(size_t index) except +yodaerr

        vector[double] xEdges() except +yodaerr
        vector[double] xWidths() except +yodaerr

        double xMin() except +yodaerr
        double xMax() except +yodaerr

        size_t numBins() except +yodaerr
        size_t numBinsX() except +yodaerr

        vector[ProfileBin1D] bins() #except +yodaerr
        int binIndexAt(double x) except +yodaerr
        const ProfileBin1D& bin(size_t ix) #except +yodaerr
        const ProfileBin1D& binAt(double x) #except +yodaerr

        # The trick here is to treat these not as references.
        # I suppose when you think about it, it makes sense
        Dbn2D& totalDbn()
        Dbn2D& underflow()
        Dbn2D& overflow()

        unsigned long numEntries(bool)
        double effNumEntries(bool)
        double sumW(bool)
        double sumW2(bool)

        double xMean(bool)
        double xVariance(bool)
        double xStdDev(bool)
        double xStdErr(bool)
        double xRMS(bool)

        operator + (Profile1D)
        operator - (Profile1D)
        operator / (Profile1D)

    Scatter2D Profile1D_div_Profile1D "divide" (const Profile1D&, const Profile1D&) except +yodaerr

cdef extern from "merge.hh":
    void Profile1D_iadd_Profile1D "cython_iadd" (Profile1D*, Profile1D*)
    void Profile1D_isub_Profile1D "cython_isub" (Profile1D*, Profile1D*)
    # void Profile1D_imul_dbl "cython_imul_dbl" (Profile1D*, double)
    # void Profile1D_idiv_dbl "cython_idiv_dbl" (Profile1D*, double)
    Profile1D* Profile1D_add_Profile1D "cython_add" (Profile1D*, Profile1D*)
    Profile1D* Profile1D_sub_Profile1D "cython_sub" (Profile1D*, Profile1D*)
    Profile1D* Profile1D_div_Profile1D "cython_div" (Profile1D*, Profile1D*)

cdef extern from "YODA/Scatter2D.h" namespace "YODA":
    Scatter2D mkScatter_Profile1D "YODA::mkScatter" (const Profile1D&, bool, bool, double, double) except +yodaerr

#}}} Profile1D



# Profile2D {{{
cdef extern from "YODA/Profile2D.h" namespace "YODA":
    cdef cppclass Profile2D(AnalysisObject,Fillable,Binned):
        Profile2D() except +yodaerr

        Profile2D(string path, string title) except +yodaerr

        Profile2D(size_t nbinsX, double lowerX, double upperX,
                  size_t nbinsY, double lowerY, double upperY,
                  string path, string title) except +yodaerr

        Profile2D(vector[double] xedges, vector[double] yedges,
                  string path, string title) except +yodaerr

        Profile2D(Profile2D p, string path) except +yodaerr

        #Profile2D(Scatter3D s, string path) except +yodaerr

        #Profile2D(Histo2D p, string path)

        Profile2D clone() except +yodaerr
        Profile2D* newclone() except +yodaerr

        # TODO: add missing functions and enable refs + exceptions when Cython allows

        void reset() except +yodaerr

        void fill(double x, double y, double z, double weight, double fraction) except +yodaerr
        void fillBin(size_t i, double z, double weight, double fraction) except +yodaerr

        void scaleW(double s) except +yodaerr
        void scaleXY(double, double)

        # void mergeBins(size_t, size_t) except +yodaerr
        # void rebin(unsigned int n) except +yodaerr

        size_t numBins() except +yodaerr
        size_t numBinsX() except +yodaerr
        size_t numBinsY() except +yodaerr

        vector[ProfileBin2D]& bins() #except +yodaerr
        int binIndexAt(double x, double y) except +yodaerr
        const ProfileBin2D& bin(size_t ix) #except +yodaerr
        const ProfileBin2D& binAt(double x, y) #except +yodaerr

        void addBin(const pair[double, double]&, const pair[double, double]&) except +yodaerr
        void addBins(const vector[double]&, const vector[double]&) except +yodaerr
        void rmBin(size_t index) except +yodaerr

        vector[double] xEdges() except +yodaerr
        vector[double] yEdges() except +yodaerr
        vector[double] xWidths() except +yodaerr
        vector[double] yWidths() except +yodaerr

        double xMin() except +yodaerr
        double xMax() except +yodaerr
        double yMin() except +yodaerr
        double yMax() except +yodaerr

        # Dbn3D& outflow(int, int) #except +yodaerr

        # Whole histo data
        Dbn3D& totalDbn() #except +yodaerr

        unsigned long numEntries(bool)
        double effNumEntries(bool)
        double sumW(bool)
        double sumW2(bool)

        double xMean(bool)
        double yMean(bool)
        double xVariance(bool)
        double yVariance(bool)
        double xStdDev(bool)
        double yStdDev(bool)
        double xStdErr(bool)
        double yStdErr(bool)
        double xRMS(bool)
        double yRMS(bool)

        operator + (Profile2D)
        operator - (Profile2D)
        operator / (Profile2D)

    Scatter3D Profile2D_div_Profile2D "divide" (const Profile2D&, const Profile2D&) except +yodaerr

cdef extern from "merge.hh":
    void Profile2D_iadd_Profile2D "cython_iadd" (Profile2D*, Profile2D*)
    void Profile2D_isub_Profile2D "cython_isub" (Profile2D*, Profile2D*)
    # void Profile2D_imul_dbl "cython_imul_dbl" (Profile2D*, double)
    # void Profile2D_idiv_dbl "cython_idiv_dbl" (Profile2D*, double)
    Profile2D* Profile2D_add_Profile2D "cython_add" (Profile2D*, Profile2D*)
    Profile2D* Profile2D_sub_Profile2D "cython_sub" (Profile2D*, Profile2D*)
    Profile2D* Profile2D_div_Profile2D "cython_div" (Profile2D*, Profile2D*)

cdef extern from "YODA/Scatter3D.h" namespace "YODA":
    Scatter3D mkScatter_Profile2D "YODA::mkScatter" (const Profile2D&, bool, bool) except +yodaerr

#}}} Profile2D





# Streams {{{

cdef extern from "<iostream>" namespace "std":
    cdef cppclass istream:
        istringstream()
        string& str(string&)

cdef extern from "<sstream>" namespace "std":
    cdef cppclass istringstream:
        istringstream()
        string& str(string&)
    cdef cppclass ostringstream:
        ostringstream()
        string& str()


cdef extern from "YODA/IO.h" namespace "YODA":
    void IO_read_from_file "YODA::read" (string&, vector[AnalysisObject*]&) except +yodaerr
    void IO_read_from_stream "YODA::read" (istream&, vector[AnalysisObject*]& aos, string&) except +yodaerr
    void IO_read_from_stringstream "YODA::read" (istringstream&, vector[AnalysisObject*]& aos, string&) except +yodaerr

cdef extern from "YODA/Index.h" namespace "YODA":
    cdef cppclass Index:
        unordered_map[string, unordered_map[string, int]] getAOIndex() except +yodaerr
        string toString() except +yodaerr

cdef extern from "YODA/Reader.h" namespace "YODA":
    cdef cppclass Reader:
        void read(istringstream&, vector[AnalysisObject*]&) except +yodaerr
        void read_from_file "YODA::Reader::read" (string&, vector[AnalysisObject*]&) except +yodaerr
        Index& make_index "YODA::Reader::mkIndex" (string&)

cdef extern from "YODA/ReaderYODA.h" namespace "YODA":
    Reader& ReaderYODA_create "YODA::ReaderYODA::create" ()

cdef extern from "YODA/ReaderFLAT.h" namespace "YODA":
    Reader& ReaderFLAT_create "YODA::ReaderFLAT::create" ()

cdef extern from "YODA/ReaderAIDA.h" namespace "YODA":
    Reader& ReaderAIDA_create "YODA::ReaderAIDA::create" ()

cdef extern from "YODA/Reader.h" namespace "YODA":
    Reader& Reader_create "YODA::mkReader" (string& filename)


cdef extern from "YODA/IO.h" namespace "YODA":
    void IO_write_to_file "YODA::write" (string&, vector[AnalysisObject*]&, int) except +yodaerr

cdef extern from "YODA/Writer.h" namespace "YODA":
    cdef cppclass Writer:
        void write(ostringstream&, vector[AnalysisObject*]&) except +yodaerr
        void write_to_file "YODA::Writer::write" (string&, vector[AnalysisObject*]&) except +yodaerr
        void setPrecision(int precision)
        void useCompression(bool compress)

cdef extern from "YODA/WriterYODA.h" namespace "YODA":
    Writer& WriterYODA_create "YODA::WriterYODA::create" ()

cdef extern from "YODA/WriterFLAT.h" namespace "YODA":
    Writer& WriterFLAT_create "YODA::WriterFLAT::create" ()

cdef extern from "YODA/WriterAIDA.h" namespace "YODA":
    Writer& WriterAIDA_create "YODA::WriterAIDA::create" ()

cdef extern from "YODA/Reader.h" namespace "YODA":
    Writer& Writer_create "YODA::mkWriter" (string& filename)

# Streams }}}



# Axis1D {{{
cdef extern from "YODA/Axis1D.h" namespace "YODA":
    cdef cppclass Axis1D[BIN1D, DBN]:
        Axis1D() except +yodaerr
        Axis1D(vector[double] binedges) except +yodaerr
        Axis1D(size_t, double, double) except +yodaerr
        Axis1D(vector[BIN1D] bins) except +yodaerr
        void addBin(double, double) except +yodaerr
        size_t numBins() except +yodaerr
        vector[BIN1D]& bins()
        double xMin() except +yodaerr
        double xMax() except +yodaerr
        vector[double] xEdges() except +yodaerr
        long getBinIndex(double)
        void reset()
        DBN& totalDbn()
        DBN& underflow()
        DBN& overflow()
        void eraseBin(size_t index) except +yodaerr
        void mergeBins(size_t, size_t) except +yodaerr
# Axis1D }}}


# Axis2D {{{
cdef extern from "YODA/Axis2D.h" namespace "YODA":
    cdef cppclass Axis2D[BIN2D, DBN]:
        Axis2D() except +yodaerr
        Axis2D(vector[double], vector[double]) except +yodaerr
        Axis2D(size_t, pair[double, double], size_t, pair[double, double]) except +yodaerr
        Axis2D(vector[BIN2D] bins) except +yodaerr
        void addBin(pair[double, double], pair[double, double]) except +yodaerr
        size_t numBins() except +yodaerr
        vector[BIN2D]& bins()
        double xMin() except +yodaerr
        double xMax() except +yodaerr
        double yMin() except +yodaerr
        double yMax() except +yodaerr
        long getBinIndex(double, double)
        void reset()
        DBN& totalDbn()
        # TODO: reinstate DBN& outflow(int, int)
        void eraseBin(size_t index) except +yodaerr
        void mergeBins(size_t, size_t) except +yodaerr
# Axis2D }}}
