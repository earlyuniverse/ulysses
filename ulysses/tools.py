import cmath
import ulysses

def readConfig(fname):
    from collections import OrderedDict

    fixed  = OrderedDict()
    ranges = OrderedDict()
    with open(fname) as f:
        for line in f:
            l=line.strip()
            if l.startswith("#"):
                continue
            if len(l) == 0:
                continue

            fields = l.split()
            if len(fields)==1:
                print( "Warning, not understood instruction:", l)
                continue

            elif len(fields)==2:
                try:
                    fixed[fields[0]] = float(fields[1])
                except:
                    fixed[fields[0]] = complex(fields[1])

            elif len(fields)==3:
                ranges[fields[0]] = (float(fields[1]), float(fields[2]))

            else:
                print ("Warning, not understood instruction:", l)
                continue
        return ranges, fixed


def getBuiltInModels():
    """
    Return as dict of short name a and class.
    """
    classes    = [cls                for cls in ulysses.ULSBase.__subclasses__()]
    shortnames = [cls.shortname(cls) for cls in ulysses.ULSBase.__subclasses__()]
    return {n : c for n, c in zip(shortnames, classes) if n!=""}

def selectModel(model, **kwargs):
    r"""

    This function loads and returns an instance of one of the built-in models.
    The kwargs are passed on to the base class.

        :Arguments:
            * *model* (``str``) --
              The model. Must be in ["1DME", "2DME", "3DME", "1BE", "2BE", "2resonant", "3DMEsct", "3DMErhtau"]

    """

    builtin = getBuiltInModels()

    import ulysses
    if not ":" in model:
        if model in builtin:
            return builtin[model](**kwargs)
        else:
            raise Exception("Specified model '{}' unknown.\n Select from: {}".format(model, sorted([k for k in builtin.keys()])))
    else:
        return loadPlugin(model, **kwargs)

def loadPlugin(model, **kwargs):
    r"""
    Plugin loader. This attemts to load and return and instance of
    the class CLASS in the file FILENAME. Both are given as a single
    string using : as separator. The kwargs are passed on to the base class.

        :Arguments:
            * *model* (``str``) --
              FILENAME:CLASS

    """
    import ulysses
    m_file, m_name = model.split(":")
    import os
    if not os.path.exists(m_file):
        raise Exception("Specified module file '{}' does not exist.{}".format(m_file))

    # https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path
    from importlib.machinery import SourceFileLoader

    foo = SourceFileLoader(m_name, m_file).load_module()

    # https://stackoverflow.com/questions/1796180/how-can-i-get-a-list-of-all-classes-within-current-module-in-python
    import inspect

    IM = dict(inspect.getmembers(foo, inspect.isclass))
    if not m_name in IM:
        raise Exception("Specified class '{}'  not found in module {}".format(m_name, m_file))
    if not issubclass(IM[m_name], ulysses.ULSBase):
        raise Exception("Specified class '{}'  does not derive from leptoms.ULSBase".format(m_name))

    return IM[m_name](**kwargs)
