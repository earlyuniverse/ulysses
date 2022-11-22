#! /usr/bin/env python

from __future__ import print_function
import yoda, random

h = yoda.Histo2D(5,0.,10., 5,0.,10., "/foo")
for _ in range(100):
    h.fill(random.gauss(5, 3), random.gauss(5, 2))
print(h)
print(h.bins())
print(h.bin(2))

yoda.write([h], "h2d.yoda")
aos = yoda.read("h2d.yoda")
for _, ao in aos.items():
    print(ao)

yoda.write([h], "h2d.dat")
# aos = yoda.read("h2d.dat")
# for _, ao in aos.iteritems():
#     print ao

# Check scatter conversion
s = yoda.mkScatter(h)
s = h.mkScatter()
s2 = s.mkScatter()

if h.numBins() != s.numPoints():
    print("FAIL mkScatter() #bin={} -> #point={}".format(h.numBins(), s.numPoints()))
    exit(11)
if h.yVals()[0] != s.point(0).y():
    print("FAIL mkScatter() bin0 value={} -> point0 value={}".format(h.yVal(0), s.point(0).y()))
    exit(12)
s = yoda.mkScatter(h, h_binsizediv=False)
if h.sumWs()[0] != s.point(0).z():
    print("FAIL mkScatter(h_binsizediv=False) bin0 value={} -> point0 value={}".format(h.sumWs()[0], s.point(0).z()))
    exit(21)


def test_addBins():
    h = yoda.Histo2D("/bar")
    NBINS_X = 100000
    j = 0
    h.addBins([(i, i + 1, j, j + 1) for i in range(NBINS_X)])
    j += 1
    h.addBins([(float(i), float(i + 1), float(j), float(j + 1)) for i in range(NBINS_X)])
    j += 1
    h.addBins([yoda.HistoBin2D(i, i + 1, j, j + 1) for i in range(NBINS_X)])
    j += 1
    h.addBins(iter([(i, i + 1, j, j + 1) for i in range(NBINS_X)]))
    j += 1
    h.addBins(iter([yoda.HistoBin2D(i, i + 1, j, j + 1) for i in range(NBINS_X)]))
    j += 1
    h.fillBin(500, 123.0)
    h.fillBin(NBINS_X * 1 + 500, 321.0)
    assert len(h.bins()) == NBINS_X * j
    assert h.totalDbn().sumW() == 123.0 + 321.0

test_addBins()
