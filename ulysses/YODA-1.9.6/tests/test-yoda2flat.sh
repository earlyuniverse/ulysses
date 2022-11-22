#!/bin/bash

set -e

yoda2flat ${YODA_TESTS_SRC}/test1.yoda yoda2flat.dat

diff -u yoda2flat.dat ${YODA_TESTS_SRC}/yoda2flat-ref.dat

rm -f yoda2flat.dat
