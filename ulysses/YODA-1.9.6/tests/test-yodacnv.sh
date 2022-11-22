#!/bin/bash

set -e

yodacnv ${YODA_TESTS_SRC}/test1.yoda yodacnv.dat

diff -u yodacnv.dat ${YODA_TESTS_SRC}/yoda2flat-ref.dat

rm -f yodacnv.dat
