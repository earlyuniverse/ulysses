#! /usr/bin/env python

import yoda, random

p = yoda.Profile2D(5,0.,10., 5,0.,10., "/bar")
for _ in range(10000):
    p.fill(random.gauss(1, 3), random.gauss(1, 2), random.gauss(10, 0.5))
print(p)

yoda.write([p], "p2d.yoda")
aos = yoda.read("p2d.yoda")
for _, ao in aos.items():
    print(ao)

yoda.write([p], "p2d.dat")
# aos = yoda.read("p2d.dat")
# for _, ao in aos.iteritems():
#     print ao

## Check scatter conversion
s = yoda.mkScatter(p)
s = p.mkScatter()
s2 = s.mkScatter()

if p.numBins() != s.numPoints():
    print("FAIL mkScatter() #bin={} -> #point={}".format(p.numBins(), s.numPoints()))
    exit(11)
if p.yVals()[0] != s.point(0).y():
    print("FAIL mkScatter() bin0 value={} -> point0 value={}".format(p.yVal(0), s.point(0).y()))
    exit(12)
if p.bin(0).stdErr() != s.point(0).zErrAvg():
    print("FAIL mkScatter(p_usestddev=False) bin0 err={} -> point0 err={}".format(p.bin(0).stdErr(), s.point(0).zErrAvg()))
    exit(13)
s = yoda.mkScatter(p, p_usestddev=True)
if p.bin(0).stdDev() != s.point(0).zErrAvg():
    print("FAIL mkScatter(p_usestddev=True) bin0 err={} -> point0 err={}".format(p.bin(0).stdDev(), s.point(0).zErrAvg()))
    exit(21)
