from __future__ import print_function

## Pull in core YODA C++/Python extension functionality
from yoda.core import *

__version__ = core.version()


## Pull in useful tools from submodules
from yoda.search import match_aos


## Pull in plotting tools (requires matplotlib and numpy)
try:
    import yoda.plotting
    from yoda.plotting import mplinit, plot
    HAS_PLOTTING = True
    # def plot(*args, **kwargs):
    #     from yoda.plotting import plot as p
    #     return p(*args, **kwargs)
except:
    HAS_PLOTTING = False


## Try to pull in optional ROOT compatibility
try:
    import yoda.root
    HAS_ROOT_SUPPORT = True
    # TODO: remove in v2
    def to_root(ao):
        print("yoda.to_root() is deprecated: use yoda.root.to_root()")
        return yoda.root.to_root(ao)
except:
    HAS_ROOT_SUPPORT = False
