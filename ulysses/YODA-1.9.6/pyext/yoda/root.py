# -*- python -*-

from __future__ import print_function

try:
    from yoda.rootcompat import *
except ImportError as e:
    #print("YODA built without ROOT support")
    raise e

def getall(d, basepath="/", verbose=False):
    "Recurse through a ROOT file/dir and generate (path, obj) pairs"
    try:
        for key in d.GetListOfKeys():
            kname = key.GetName()
            if key.IsFolder():
                # TODO: -> "yield from" in Py3
                for i in getall(d.Get(kname), basepath+kname+"/"):
                    yield i
            else:
                yield basepath+kname, d.Get(kname)
    except AttributeError:
        # deal with things like TObjArray that don't have GetListOfKeys()
        dname  = 'unknownName'
        dtitle = 'unknownTitle'
        dclass = 'unknownClass'
        try:
            dname  = d.GetName()
            dtitle = d.GetTitle()
            try:
                dclass = d.ClassName()
            except:
                dclass = "unknownClassType"
        except:
            dclass = 'failed GetName/GetTitle'
        if ( verbose ):
            sbase='failed d.GetListOfKeys() for name={} title="{}" class={}'
            print(sbase.format(dname,dtitle,dclass))
