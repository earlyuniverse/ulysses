#! /usr/bin/env python

import re


## A tool for filtering AO collections by path patterns
def match_aos(aos, patts, unpatts=None, search=False):
    """
    Filter a list of analysis objects to those which match given path-matching patterns.

    @a patts is a regex or iterable of regexes for positive matching, i.e. retention;
    @a unpatts is the equivalent for negative matching, i.e. rejection even if a patt matches.

    @a search will use Python regex search mode rather than match mode, i.e. match if any
    part of the path fits the regex, rather than requiring a match from the start of the path.
    """

    rtn = []

    ## Internally we just use the list of AOs, not dict keys
    if type(aos) is dict:
        aos = list(aos.values()) #< TODO: maybe iter would be ok?

    ## Normalise the pattern arg inputs
    if patts and type(patts) is str:
        patts = [patts]
    re_patts = [re.compile(p) for p in patts] if patts else None
    if unpatts and type(unpatts) is str:
        unpatts = [unpatts]
    re_unpatts = [re.compile(p) for p in unpatts] if unpatts else None

    ## Apply pattern matching to each AO
    for ao in aos:
        match = False
        if re_patts:
            for rp in re_patts:
                if (not search and rp.match(ao.path())) or (search and rp.search(ao.path())):
                    match = True
                    break
        if match and re_unpatts:
            for rp in re_unpatts:
                if (not search and rp.match(ao.path())) or (search and rp.search(ao.path())):
                    match = False
                    break
        if match:
            rtn.append(ao)

    return rtn



class PointMatcher(object):
    """\
    System for selecting subsets of bins based on a search range
    syntax extended from Professor weight files:

    Path structure: /path/parts/to/histo[syst_variation]@xmin:xmax
                or: /path/parts/to/histo[syst_variation]#nmin:nmax

    TODO: Extend to multi-dimensional ranges i.e. @xmin:xmax,#nymin:nymax,...
    """

    def __init__(self, patt):
        self.re_patt = re.compile(r"([^#@]+)(#\d+|@[\d\.:]+)?")
        self.set_patt(patt)

    def set_patt(self, patt):
        "Find path and index/pos parts of patt and assign them to object attrs"
        self.patt = None
        self.path = None
        self.indextype = None
        self.index = None
        if not patt:
            return
        ## Strip separated comments
        patt = re.sub(r"(^|\s+)#.*", "", patt)
        self.patt = patt.strip()
        match = self.re_patt.match(self.patt)
        if match:
            self.path = re.compile(match.group(1))
            if match.group(2):
                # TODO: handle mixed-type ranges?
                self.indextype, indexstr = match.group(2)[0], match.group(2)[1:]
                if self.indextype:
                    if not ":" in indexstr:
                        self.index = float(indexstr)
                    else:
                        indexstr2 = indexstr.split(":", 1)
                        if not indexstr2[0]: indexstr2[0] = "-inf"
                        if not indexstr2[1]: indexstr2[1] = "inf"
                        self.index = [float(istr) for istr in indexstr2]

    def match_path(self, path):
        return self.path.match(path) is not None

    def search_path(self, path):
        return self.path.search(path) is not None

    def match_pos(self, p):
        """Decide if a given point p is in the match range.

        p must be an object with attrs xmin, xmax, n

        TODO: Use open ranges to include underflow and overflow

        TODO: Allow negative indices in Python style, and use index=-1
        to mean the N+1 index needed to include the last bin without
        picking up the overflow, too.

        TODO: Extension to multiple dimensions
        """
        if not self.indextype:
            accept = True
        elif self.indextype == "#":
            if type(self.index) is float:
                accept = (p.n == int(self.index))
            else:
                accept = (p.n >= self.index[0] and p.n < self.index[1])
        else: # self.indextype == "@"
            if type(self.index) is float:
                accept = (self.index >= p.xmin and self.index < p.xmax)
            else:
                accept = (p.xmax > self.index[0] and p.xmin <= self.index[1])
        return accept

    def __repr__(self):
        s = "PointMatcher('%s' %s %s %s)" % (self.patt, self.path, self.indextype, self.index)
        return s


if __name__ == "__main__":

    from . import linspace, logspace
    from pprint import pprint
    import math

    class Point(object):
        def __init__(self, path, n, xmin, xmax, value=None):
            self.path = path
            self.n = n #< bin index
            self.xmin = xmin
            self.xmax = xmax
            self.value = value
        def __repr__(self):
            val = " = {}".format(self.value) if self.value else ""
            return "Point({} #{:d}  {:.2e}--{:.2e}{})".format(self.path, self.n, self.xmin, self.xmax, val)

    POINTS = []
    xs = linspace(5, 0, 10)
    POINTS += [Point("/foo", i, xs[i], xs[i+1]) for i in range(len(xs)-1)]
    xs = linspace(50, 0, 100)
    POINTS += [Point("/bar", i, xs[i], xs[i+1]) for i in range(len(xs)-1)]
    xs = linspace(20, 0, math.pi)
    POINTS += [Point("/baz/pi", i, xs[i], xs[i+1]) for i in range(len(xs)-1)]
    xs = logspace(20, 0.1, 50)
    POINTS += [Point("/baz/log", i, xs[i], xs[i+1]) for i in range(len(xs)-1)]
    pprint(POINTS)

    pms = [PointMatcher("/foo   # this bit is a comment"),
           PointMatcher("/bar#1  # this bit is a comment"),
           PointMatcher("/baz/pi@2.76  # this bit is a comment"),
           PointMatcher("/baz/.*@:1.32 # this bit is a comment"),
           PointMatcher("/baz/.*@2.76: # this bit is a comment")]
    pprint(pms)

    for pm in pms:
        print()
        print(pm.patt, pm.indextype, pm.index)
        pprint( [p for p in POINTS if pm.match_path(p.path) and pm.match_pos(p)] )
