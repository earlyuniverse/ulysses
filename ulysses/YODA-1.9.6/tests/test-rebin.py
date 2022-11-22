#! /usr/bin/env python

import yoda
import numpy as np

h = yoda.Histo1D(10, 0, 5)
for x in (0.1, 0.2, 1.3, 2.1, 2.7, 2.8, 4.0, 4.1):
    h.fill(x)
print(h.numBins())
print(h.xEdges())

print("Ha1")
ha1 = h.clone()
ha1.rebinBy(2)
print(ha1.numBins)
print(ha1.xEdges())
assert ha1.numBins() == 5
assert np.allclose(ha1.xEdges(), [0.0, 1.0, 2.0, 3.0, 4.0, 5.0])

print("Ha2")
ha2 = h.clone()
ha2.rebinBy(2, 2, 7)
print(ha2.numBins())
print(ha2.xEdges())
assert ha2.numBins() == 7
assert np.allclose(ha2.xEdges(), [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 4.5, 5.0])

print("Hb1")
hb1 = h.clone()
hb1.rebinTo([0., 1., 3., 5.])
print(hb1.numBins())
print(hb1.xEdges())
assert hb1.numBins() == 3
assert np.allclose(hb1.xEdges(), [0.0, 1.0, 3.0, 5.0])

print("Hb2")
hb2 = h.clone()
hb2.rebin([1., 1.5,  3., 4.5])
print(hb2.numBins())
print(hb2.xEdges())
assert hb2.numBins() == 3
assert np.allclose(hb2.xEdges(), [1.0, 1.5, 3.0, 4.5])
