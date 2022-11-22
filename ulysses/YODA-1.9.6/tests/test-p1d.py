#! /usr/bin/env python

import yoda, random

p1 = yoda.Profile1D(20, 0.0, 100.0, path="/foo", title="MyTitle")

linspace = yoda.linspace(20, 0.0, 100.0)
p2 = yoda.Profile1D(linspace, path="/bar", title="Linearly spaced histo")

logspace = yoda.logspace(20, 1.0, 64)
p3 = yoda.Profile1D(logspace, path="/baz", title="Log-spaced histo")


NUM_SAMPLES = 1000
for i in range(NUM_SAMPLES):
    # exp = - (i-NUM_SAMPLES/2)**2 / float(NUM_SAMPLES/4)
    # val = 2.718 ** exp
    val = random.uniform(0,100)
    p1.fill(val, random.gauss(5, 3));
    p2.fill(val, random.gauss(5, 4));
    p3.fill(val, random.gauss(5, 5));
print(p1)
print(p2)
print(p3)


yoda.write([p1,p2,p3], "p1d.yoda")
aos = yoda.read("p1d.yoda")
for _, ao in aos.items():
    print(ao)

yoda.writeFLAT([p1,p2,p3], "p1d.dat")
aos = yoda.read("p1d.dat")
for _, ao in aos.items():
    print(ao)

s = yoda.mkScatter(p1)
s = p1.mkScatter()
s2 = s.mkScatter()


# Check that the bin scaling is done properly
s1 = yoda.mkScatter(p1)
if p1.numBins() != s1.numPoints():
    print("FAIL mkScatter() #bin={} -> #point={}".format(p1.numBins(), s1.numPoints()))
    exit(11)
if p1.yVals()[0] != s1.point(0).y():
    print("FAIL mkScatter() bin0 value={} -> bin0 value={}".format(p1.yVal(0), s1.point(0).y()))
    exit(12)
if p1.yErrs(sd=False)[0] != s2.point(0).errAvg(2):
    print("FAIL mkScatter() bin0 err={} -> bin0 err={}".format(p1.yErrs(sd=False)[0], s2.errAvg(0).y()))
    exit(22)

# Check that the bin scaling is done properly
s2 = yoda.mkScatter(p1, p_usestddev=True)
if p1.numBins() != s2.numPoints():
    print("FAIL mkScatter(h_binsizediv=True) #bin={} -> #point={}".format(p1.numBins(), s1.numPoints()))
    exit(21)
if p1.yErrs(sd=True)[0] != s2.point(0).errAvg(2):
    print("FAIL mkScatter(h_binsizediv=True) bin0 err={} -> point0 err={}".format(p1.yErrs(sd=True)[0], s2.point(0).yErrAvg()))
    exit(22)

## Check the inclusion of underflow and overflow bins
su = yoda.mkScatter(p1, uflow_binwidth=1.0)
so = yoda.mkScatter(p1, oflow_binwidth=1.0)
suo = yoda.mkScatter(p1, uflow_binwidth=1.0, oflow_binwidth=1.0)
if p1.numBins() != (su.numPoints()-1):
    print("FAIL mkScatter(uflow) #bin={} -> #point={}".format(p1.numBins(), su.numPoints()))
    exit(31)
if p1.numBins() != (so.numPoints()-1):
    print("FAIL mkScatter(oflow) #bin={} -> #point={}".format(p1.numBins(), so.numPoints()))
    exit(32)
if p1.numBins() != (suo.numPoints()-2):
    print("FAIL mkScatter(uflow, oflow) #bin={} -> #point={}".format(p1.numBins(), suo.numPoints()))
    exit(33)
