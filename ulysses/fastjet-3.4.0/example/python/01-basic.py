#!/usr/bin/env python
"""Simple example to try out fastjet from python, with funcionality similar to
../01-basic.cc

For this script to work, make sure that the installation location for
the fastjet python module (cf. fastjet-config --pythonpath) is
included in your PYTHONPATH environment variable.

"""
from __future__ import print_function

import fastjet as fj
#import gzip

def main():


    # set up our jet definition and a jet selector
    jet_def = fj.JetDefinition(fj.antikt_algorithm, 0.6)
    selector = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(4.5)

    filename = '../data/single-event.dat'
    
    # get the event
    event = read_event(filename)
    # cluster it
    jets = selector(jet_def(event))

    # print out some information about the event and clustering
    print("Event has {0} particles".format(len(event)))
    print("jet definition is:",jet_def)
    print("jet selector is:", selector,"\n")
    
    # print the jets
    print_jets(jets)
    
    # get internal information about one of the jets
    if (len(jets) > 0):
        print("Number of constituents of jets[0] is {0}".format(len(jets[0].constituents())))

#----------------------------------------------------------------------
def print_jets(jets):
    print("{0:>5s} {1:>10s} {2:>10s} {3:>10s}".format("jet #", "pt", "rap", "phi"))
    for ijet in range(len(jets)):
        print("{0:5d} {1:10.3f} {2:10.4f} {3:10.4f}".format(
            ijet, jets[ijet].pt(), jets[ijet].rap(), jets[ijet].phi()))
    
        
#----------------------------------------------------------------------
def read_event(filename):
    f = open(filename, 'r')
    event = []
    for line in f:
        p = line.split()
        event.append(fj.PseudoJet(float(p[0]),float(p[1]),float(p[2]),float(p[3])));

    return event
    
if __name__ == '__main__':
    main()
