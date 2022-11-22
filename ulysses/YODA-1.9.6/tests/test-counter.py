#! /usr/bin/env python

import yoda, random

c = yoda.Counter(path="/foo", title="MyTitle")

NUM_SAMPLES = 1000
for i in range(NUM_SAMPLES):
    c.fill(random.gauss(10,3))
print(c.val(), "+-", c.err())

yoda.write([c], "counter.yoda")
aos = yoda.read("counter.yoda")
for _, ao in aos.items():
    print(ao)
