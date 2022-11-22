#! /usr/bin/env python

import yoda, random

h1 = yoda.Histo1D(20, 0.0, 100.0, path="/foo", title="MyTitle")

linspace = yoda.linspace(20, 0.0, 100.0)
h2 = yoda.Histo1D(linspace, path="/bar", title="Linearly spaced histo")

logspace = yoda.logspace(20, 1.0, 64)
h3 = yoda.Histo1D(logspace, path="/baz", title="Log-spaced histo")


NUM_SAMPLES = 1000
for i in range(NUM_SAMPLES):
    exp = - (i-NUM_SAMPLES/2)**2 / float(NUM_SAMPLES/4)
    val = 2.718 ** exp
    h1.fill(val);
    h2.fill(val);
    h3.fill(val);
print(h1.xMean(), "+-", h1.xStdDev())
print(h1)
print(h1.bins())
print(h1.bin(2))
print(h2)
print(h3)


yoda.write([h1,h2,h3], "h1d.yoda")
aos = yoda.read("h1d.yoda")
for _, ao in aos.items():
    print(ao)

yoda.writeFLAT([h1,h2,h3], "h1d.dat")
aos = yoda.read("h1d.dat")
for _, ao in aos.items():
    print(ao)
s = yoda.mkScatter(h1)
s = h1.mkScatter()
s2 = s.mkScatter()


# Check that the bin scaling is done properly
s1 = yoda.mkScatter(h3)
if h3.numBins() != s1.numPoints():
    print("FAIL mkScatter() #bin={} -> #point={}".format(h3.numBins(), s1.numPoints()))
    exit(11)
if h3.yVals()[0] != s1.point(0).y():
    print("FAIL mkScatter() bin0 value={} -> bin0 value={}".format(h3.yVal(0), s1.point(0).y()))
    exit(12)

# Check that the bin scaling is done properly
s2 = yoda.mkScatter(h3, h_binsizediv=False)
if h3.numBins() != s2.numPoints():
    print("FAIL mkScatter(h_binsizediv=False) #bin={} -> #point={}".format(h3.numBins(), s1.numPoints()))
    exit(21)
if h3.yVals(area=True)[0] != s2.point(0).val(2):
    print("FAIL mkScatter(h_binsizediv=False) bin0 value={} -> point0 value={}".format(h3.yVals(area=True)[0], s2.point(0).y()))
    exit(22)

## Check the inclusion of underflow and overflow bins
su = yoda.mkScatter(h3, uflow_binwidth=1.0)
so = yoda.mkScatter(h3, oflow_binwidth=1.0)
suo = yoda.mkScatter(h3, uflow_binwidth=1.0, oflow_binwidth=1.0)
if h3.numBins() != (su.numPoints()-1):
    print("FAIL mkScatter(uflow) #bin={} -> #point={}".format(h3.numBins(), su.numPoints()))
    exit(31)
if h3.numBins() != (so.numPoints()-1):
    print("FAIL mkScatter(oflow) #bin={} -> #point={}".format(h3.numBins(), so.numPoints()))
    exit(32)
if h3.numBins() != (suo.numPoints()-2):
    print("FAIL mkScatter(uflow, oflow) #bin={} -> #point={}".format(h3.numBins(), suo.numPoints()))
    exit(33)
