#! /usr/bin/env python

import yoda

hs = [yoda.Histo1D(5, 0, 10, "/H%d" % i) for i in range(3)]

hs[0].fill(2)
print(hs[0])
print()

hs[1].fill(5)
print(hs[1])
hs[1].scaleW(2)
print(hs[1])
print()

hs[2].fill(8, 2)
print(hs[2])
hs[2] += hs[1]
print(hs[2])
print()

print(hs[0].annotations)
print(hs[1].annotations)

h = hs[0] + hs[1]
print(hs[0])
print(hs[1])
print(h)

# TODO: Not currently supported... some Cython problem with type overloading in __add__, __radd__
# h = sum(hs)
# print h
