#!/usr/bin/env python
"""Simple example to try out fastjet from python, showing multi-event reading.

For this script to work, make sure that the installation location for
the fastjet python module (cf. fastjet-config --pythonpath) is
included in your PYTHONPATH environment variable.

"""
from __future__ import print_function

import fastjet as fj
#import gzip

def main():

    # get the banner out of the way early on
    fj.ClusterSequence.print_banner()
    print()

    # set up our jet definition and a jet selector
    jet_def = fj.JetDefinition(fj.antikt_algorithm, 0.4)
    selector = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(4.5)
    print("jet definition is:",jet_def)
    print("jet selector is:", selector,"\n")

    #filename = '../data/single-event.dat'
    filename = '../data/Pythia-PtMin1000-LHC-10ev.dat'
    f = open(filename,'r')
    #filename = '/Users/gsalam/work/fastjet/data/Pythia-PtMin50-LHC-10kev.dat.gz'
    #f = gzip.GzipFile(filename,'rb')
    
    # get the event
    iev = 0
    while True:
        event = read_event(f)
        iev += 1
        if (len(event) == 0): break
        jets = selector(jet_def(event))
        print("Event {0} has {1} particles".format(iev, len(event)))
        
        # cluster it
        for ijet in range(len(jets)):
            print("jet {0} pt and rap: {1} {2}".format(ijet, jets[ijet].pt(), jets[ijet].rap()))
            
        # make sure jet-related information is correctly held
        if (len(jets) > 0):
            print("Number of constituents of jets[0] is {0}".format(len(jets[0].constituents())))
            
#----------------------------------------------------------------------
def read_event(file_or_filename):
    """
Routine that can take either an existing opened file object, or a
filename (which it will open itself) and then reads an event from that
file. An event is deemed to end when the file comes to an end or when
the reader encounters the string "#END".

The event is converted to a python list of PseudoJets
    """

    # open the file if necessary
    if (isinstance(file_or_filename,str)) : f = open(file_or_filename, 'r')
    else                                  : f = file_or_filename

    # create an empty list
    event = []
    while True:
        line = f.readline()
        # exit if the file has come to an end
        if   (not line): break
        # or if we reach the string "#END"
        if   (len(line) >=4 and line[0:4] == '#END'): break

        # ignore comment lines or empty lines
        elif (line[0] == '#' or len(line) <= 1): continue

        # assume we have a good line and split it into px, py, pz, E
        p = line.split()
        # and append the PseudoJet
        event.append(fj.PseudoJet(float(p[0]),float(p[1]),float(p[2]),float(p[3])));

    return event
    
if __name__ == '__main__':
    main()
