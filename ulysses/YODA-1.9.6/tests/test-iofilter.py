#! /usr/bin/env python

import yoda
import os

srcpath = os.getenv('YODA_TESTS_SRC')

testfile = os.path.join(srcpath,"iofilter.yoda")

aos = yoda.read(testfile, False, patterns=r".*_(eta|mass)_(2|4)")
print(aos)
assert len(aos) == 4

index = yoda.mkIndexYODA(testfile)
print("Index of ref file:")
print(index)
import re
aos = yoda.read(testfile, False, patterns=[r".*_y_(1|3)", re.compile(r".*_dphi_(2|4)")])
print(aos)
assert len(aos) == 3

aos = yoda.read(testfile, False, patterns=[r".*_y_(1|3)", re.compile(r".*_dphi_(2|4)")], unpatterns=r".*_y_1")
print(aos)
assert len(aos) == 2
