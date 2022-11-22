#! /usr/bin/env python

import yoda, random

h1, h2 = [yoda.Histo1D(4, 0, 10) for _ in range(2)]
for i in range(1000):
    h1.fill(random.uniform(0,10))
    h2.fill(random.uniform(0,10))

s = h1 / h2
print(s)
for p in s.points():
    print(" ", p)

print()

s = h1.divideBy(h2)
print(s)
for p in s.points():
    print(" ", p)
