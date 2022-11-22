#!/bin/bash

set -e

cat ${YODA_TESTS_SRC}/yodahist-fills.txt | yodahist h1 10 0. 100. out yodahist.yoda title foobar path /foo/bar/baz
yodadiff yodahist.yoda ${YODA_TESTS_SRC}/yodahist-ref.yoda

rm -f yodahist.yoda
