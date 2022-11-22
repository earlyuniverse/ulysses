#----------------------------------------------------------------------
# $Id: read_lund_json.py 1289 2021-11-09 11:53:53Z scyboz $
#
# Copyright (c) 2018-, Frederic A. Dreyer, Keith Hamilton, Alexander Karlberg,
# Gavin P. Salam, Ludovic Scyboz, Gregory Soyez, Rob Verheyen
#
#----------------------------------------------------------------------
# This file is part of FastJet contrib.
#
# It is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
#
# It is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this code. If not, see <http://www.gnu.org/licenses/>.
#----------------------------------------------------------------------

from abc import ABC, abstractmethod
import numpy as np
from math import log, ceil, cos, sin
import json, gzip, sys

#======================================================================
class Reader(object):
    """
    Reader for files consisting of a sequence of json objects, one per
    line Any pure string object is considered to be part of a header
    (even if it appears at the end!).

    Once you have created a Reader, it is a standard Python iterable
    object. Each iteration gives the json entry for one jet's
    declustering. 
    """
    #----------------------------------------------------------------------
    def __init__(self, infile, nmax = -1):
        self.infile = infile
        self.nmax = nmax
        self.reset()

    #----------------------------------------------------------------------
    def set_nmax(self, nmax):
        self.nmax = nmax

    #----------------------------------------------------------------------
    def reset(self):
        """
        Reset the reader to the start of the file, clear the header and event count.
        """
        if self.infile.endswith('.gz'):
            self.stream = gzip.open(self.infile,'r')
        else:
            self.stream = open(self.infile,'r')
        self.n = 0
        self.header = []
        
    #----------------------------------------------------------------------
    def __iter__(self):
        # needed for iteration to work 
        return self
        
    #----------------------------------------------------------------------
    def __next__(self):
        ev = self.next_event()
        if (ev is None): raise StopIteration
        else           : return ev

    #----------------------------------------------------------------------
    def next(self): return self.__next__()
        
    #----------------------------------------------------------------------
    def next_event(self):

        # we have hit the maximum number of events
        if (self.n == self.nmax):
            #print ("# Exiting after having read nmax jet declusterings")
            return None
        
        try:
            line = self.stream.readline()
            if self.infile.endswith('.gz'):
                j = json.loads(line.decode('utf-8'))
            else:
                j = json.loads(line)
        except IOError:
            print("# got to end with IOError (maybe gzip structure broken?) around event", self.n, file=sys.stderr)
            return None
        except EOFError:
            print("# got to end with EOFError (maybe gzip structure broken?) around event", self.n, file=sys.stderr)
            return None
        except ValueError:
            print("# got to end with ValueError (empty json entry?) around event", self.n, file=sys.stderr)
            return None

        # skip this
        if (type(j) is str):
            self.header.append(j)
            return self.next_event()
        self.n += 1
        return j

#======================================================================
class Image(ABC):
    """
    Abstract base class for batch processing Lund declusterings
    """
    def __init__(self, infile, nmax=-1):
        if (type(infile) is Reader): 
            self.reader = infile
            self.reader.set_nmax(nmax)
        else:                        
            self.reader = Reader(infile, nmax)

    #----------------------------------------------------------------------
    @abstractmethod
    def process(self, event):
        pass
    
    #----------------------------------------------------------------------
    def values(self):
        """Return the values of each event in the reader."""
        res = []
        while True:
            event = self.reader.next_event()
            if event!=None:
                res.append(self.process(event))
            else:
                break
        self.reader.reset()
        return res

        
#======================================================================
class LundImage(Image):
    """
    Class to take input file (or a reader) of jet declusterings in json
    format, one json entry per line, and transform it into lund images.

    - infile: a filename or a reader
    - nmax: the maximum number of jets to process
    - npxl: the number of bins (pixels) in each dimension
    - xval: the range of x (ln 1/Delta) values
    - yval: the range of y (ln kt) values

    Once you've created the class, call values() (inherited the abstract
    base class) and you will get a python list of images (formatted as
    2d numpy arrays).
    """
    #----------------------------------------------------------------------
    def __init__(self, infile, nmax, npxl, xval = [0.0, 7.0], yval = [-3.0, 7.0]):
        Image.__init__(self, infile, nmax)
        self.npxl = npxl
        self.xmin = xval[0]
        self.ymin = yval[0]
        self.x_pxl_wdth = (xval[1] - xval[0])/npxl
        self.y_pxl_wdth = (yval[1] - yval[0])/npxl

    #----------------------------------------------------------------------
    def process(self, event):
        """Process an event and return an image of the primary Lund plane."""
        res = np.zeros((self.npxl,self.npxl))

        for declust in event:
            x = log(1.0/declust['Delta'])
            y = log(declust['kt'])
            psi = declust['psi']
            xind = ceil((x - self.xmin)/self.x_pxl_wdth - 1.0)
            yind = ceil((y - self.ymin)/self.y_pxl_wdth - 1.0)
            # print((x - self.xmin)/self.x_pxl_wdth,xind,
            #       (y - self.ymin)/self.y_pxl_wdth,yind,':',
            #       declust['delta_R'],declust['pt2'])
            if (max(xind,yind) < self.npxl and min(xind, yind) >= 0):
                res[xind,yind] += 1
        return res
    
#======================================================================
class LundDense(Image):
    """
    Class to take input file (or a reader) of jet declusterings in json
    format, one json entry per line, and reduces them to the minimal
    information needed as an input to LSTM or dense network learning.

    - infile: a filename or a reader
    - nmax: the maximum number of jets to process
    - nlen: the size of the output array (for each jet), zero padded
            if the declustering sequence is shorter; if the declustering
            sequence is longer, entries beyond nlen are discarded.

    Calling values() returns a python list, each entry of which is a
    numpy array of dimension (nlen,2). values()[i,:]  =
    (log(1/Delta),log(kt)) for declustering i.
    """
    #----------------------------------------------------------------------
    def __init__(self,infile, nmax, nlen):
        Image.__init__(self, infile, nmax)
        self.nlen      = nlen
        
    #----------------------------------------------------------------------
    def process(self, event):
        """Process an event and return an array of declusterings."""
        res = np.zeros((self.nlen, 2))
        # go over the declusterings and fill the res array
        # with the Lund coordinates
        for i in range(self.nlen):
            if (i >= len(event)):
                break
            res[i,:] = self.fill_declust(event[i])
            
        return res

    #----------------------------------------------------------------------
    def fill_declust(self,declust):
        """ 
        Create an array of size two and fill it with the Lund coordinates
        (log(1/Delta),log(kt)).  
        """
        res = np.zeros(2)
        res[0] = log(1.0/declust['Delta'])
        res[1] = log(declust['kt'])
        return res
