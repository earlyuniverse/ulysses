#cython: embedsignature=True

## System / stdlib imports
from cython.operator cimport dereference as deref, preincrement as preinc
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.map cimport map

## YODA mapping imports
cimport yoda.declarations as c
cimport yoda.util as cutil
import yoda.util as util

## Pure Python imports
from itertools import repeat
from operator import attrgetter

def version():
    "Return YODA library version as a string"
    return c.version()

include "include/Errors.pyx"
include "include/Dbn0D.pyx"
include "include/Dbn1D.pyx"
include "include/Dbn2D.pyx"
include "include/Dbn3D.pyx"
include "include/Counter.pyx"
include "include/Bin.pyx"
include "include/Bin1D.pxi"
include "include/HistoBin1D.pyx"
include "include/ProfileBin1D.pyx"
include "include/AnalysisObject.pyx"
include "include/Histo1D.pyx"
include "include/Profile1D.pyx"
include "include/Bin2D.pxi"
include "include/ProfileBin2D.pyx"
include "include/HistoBin2D.pyx"
include "include/Histo2D.pyx"
include "include/Profile2D.pyx"
include "include/Point.pyx"
include "include/Point1D.pyx"
include "include/Point2D.pyx"
include "include/Point3D.pyx"
include "include/Scatter1D.pyx"
include "include/Scatter2D.pyx"
include "include/Scatter3D.pyx"
include "include/Functions.pyx"
include "include/IO.pyx"
