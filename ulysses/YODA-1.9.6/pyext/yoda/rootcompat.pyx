cimport rootcompat as croot
import yoda
cimport yoda.declarations as cyoda
cimport yoda.util as cutil
import ROOT
cimport cython.operator.dereference as deref
from cpython.ref cimport PyObject
#ROOT.PyConfig.IgnoreCommandLineOptions = True



cdef croot.TObject* py_to_root(object pyrootobj):
    cdef PyObject* ptr = <PyObject*>pyrootobj
    return croot.py_owned_to_root(ptr)

cdef object root_to_py(croot.TObject* tobj):
    return <object> croot.root_to_py_owned(tobj)



cdef _TH1toS2(croot.TH1D* th1d, widthscale):
    return cutil.new_owned_cls(yoda.Scatter2D, croot.toNewScatter2D(th1d, widthscale))

cdef _TP1toS2(croot.TProfile* tp1):
    return cutil.new_owned_cls(yoda.Scatter2D, croot.toNewScatter2D(tp1))

# cdef _TG1toS2(croot.TGraph* tg1):
#     return cutil.new_owned_cls(yoda.Scatter2D, croot.toNewScatter2D(tg1))

cdef _TH2toS3(croot.TH2D* th2, areascale):
    return cutil.new_owned_cls(yoda.Scatter3D, croot.toNewScatter3D(th2, areascale))

# cdef _TP2toS3(croot.TProfile* tp2):
#     return cutil.new_owned_cls(yoda.Scatter3D, croot.toNewScatter3D(tp2))

# cdef _TG2toS3(croot.TGraph2D* tg2):
#     return cutil.new_owned_cls(yoda.Scatter3D, croot.toNewScatter3D(tg2))


def to_yoda(root_obj, widthscale=False):
    cdef croot.TObject* ptr = py_to_root(root_obj)
    if isinstance(root_obj, ROOT.TProfile):
        return _TP1toS2(<croot.TProfile*>ptr)
    elif isinstance(root_obj, ROOT.TH1D):
        return _TH1toS2(<croot.TH1D*>ptr, widthscale)
    # elif isinstance(root_obj, ROOT.TGraph):
    #     return _TG1toS2(<croot.TGraph*>ptr)
    elif isinstance(root_obj, ROOT.TH2D):
        return _TH2toS3(<croot.TH2D*>ptr, widthscale)
    # elif isinstance(root_obj, ROOT.TProfile2D):
    #     return _TP2toS3(<croot.TProfile2D*>ptr, widthscale)
    # elif isinstance(root_obj, ROOT.TGraph2D):
    #     return _TG2toS3(<croot.TGraph2D*>ptr)




cdef _H1toTH1D(cyoda.Histo1D* h1d, widthscale):
    return ROOT.TH1D(root_to_py(new croot.TH1D(croot.toTH1D(deref(h1d), widthscale))))

cdef _P1toTProfile(cyoda.Profile1D* p1d):
    return ROOT.TProfile(root_to_py(new croot.TProfile(croot.toTProfile(deref(p1d)))))

cdef _H2toTH2D(cyoda.Histo2D* h2d, areascale):
    return ROOT.TH2D(root_to_py(new croot.TH2D(croot.toTH2D(deref(h2d), areascale))))

# cdef _P2toTProfile2D(cyoda.Profile2D* p2d):
#     return ROOT.TProfile2D(root_to_py(new croot.TProfile2D(croot.toTProfile2D(deref(p2d)))))


cdef _S2toTGraph(cyoda.Scatter2D* s2d):
    return ROOT.TGraphAsymmErrors(root_to_py(new croot.TGraphAsymmErrors(croot.toTGraph(deref(s2d)))))

cdef _H1toTGraph(cyoda.Histo1D* h1d, usefocus, widthscale):
    return ROOT.TGraphAsymmErrors(root_to_py(new croot.TGraphAsymmErrors(croot.toTGraph(deref(h1d), usefocus, widthscale))))

cdef _P1toTGraph(cyoda.Profile1D* p1d, usefocus):
    return ROOT.TGraphAsymmErrors(root_to_py(new croot.TGraphAsymmErrors(croot.toTGraph(deref(p1d), usefocus))))


# cdef _S3toTGraph(cyoda.Scatter3D* s3d):
#     return ROOT.TGraph2D(root_to_py(new croot.TGraph2D(croot.toTGraph(deref(s3d)))))

# cdef _H2toTGraph(cyoda.Histo2D* h2d):
#     return ROOT.TGraph2D(root_to_py(new croot.TGraph2D(croot.toTGraph(deref(h2d)))))

# cdef _P2toTGraph(cyoda.Profile2D* p2d):
#     return ROOT.TGraph2D(root_to_py(new croot.TGraph2D(croot.toTGraph(deref(p2d)))))


def to_root(cutil.Base yoda_obj, asgraph=False, usefocus=False, widthscale=False):
    cdef void* ptr = yoda_obj.ptr()
    if isinstance(yoda_obj, yoda.Histo1D):
        if asgraph:
            return _H1toTGraph(<cyoda.Histo1D*> ptr, usefocus, widthscale)
        else:
            return _H1toTH1D(<cyoda.Histo1D*> ptr, widthscale)
    elif isinstance(yoda_obj, yoda.Profile1D):
        if asgraph:
            return _P1toTGraph(<cyoda.Profile1D*> ptr, usefocus)
        else:
            return _P1toTProfile(<cyoda.Profile1D*> ptr)
    elif isinstance(yoda_obj, yoda.Histo2D):
        #return _H2toTGraph(<cyoda.Histo2D*> ptr) if asgraph else _H2toTH2D(<cyoda.Histo2D*> ptr)
        return _H2toTH2D(<cyoda.Histo2D*> ptr, widthscale)
    # elif isinstance(yoda_obj, yoda.Profile2D):
    #     return _P2toTGraph(<cyoda.Profile2D*> ptr) if asgraph else _P2toTProfile2D(<cyoda.Profile2D*> ptr)
    elif isinstance(yoda_obj, yoda.Scatter2D):
        return _S2toTGraph(<cyoda.Scatter2D*> ptr)
    # elif isinstance(yoda_obj, yoda.Scatter3D):
    #     return _S3toTGraph(<cyoda.Scatter3D*> ptr)
