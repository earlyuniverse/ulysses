#! /usr/bin/env python

from io import BytesIO, StringIO

import sys
if sys.version_info <= (3,0):
    # workaround for:
    #   TypeError: initial_value must be unicode or None, not str
    from StringIO import StringIO

import os
TESTSRCDIR = os.environ.get("YODA_TESTS_SRC", ".")
def testsrcpath(fname):
    return os.path.join(TESTSRCDIR, fname)

import yoda
aos_ref = yoda.read(testsrcpath("test.yoda"))

assert len(aos_ref.keys()) > 0
print(aos_ref.keys())

ypath = testsrcpath("test.yoda")
yzpath = ypath + ".gz"

for aos in [
        yoda.read(open(ypath, "r")),
        yoda.read(StringIO(open(ypath, "r").read())),
        yoda.readYODA(ypath),
        yoda.readYODA(open(ypath, "r")),
        yoda.readYODA(StringIO(open(ypath, "r").read())),

        # Compressed
        yoda.read(yzpath),
        yoda.read(open(yzpath, "rb")),
        yoda.read(BytesIO(open(yzpath, "rb").read())),
        yoda.readYODA(yzpath),
        yoda.readYODA(open(yzpath, "rb")),
        yoda.readYODA(BytesIO(open(yzpath, "rb").read())),
]:
    print(aos.keys())
    assert set(aos.keys()) == set(aos_ref.keys())
