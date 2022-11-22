#!/bin/bash

set -e

yodals -q ${YODA_TESTS_SRC}/rivetexample.yoda -m DELPHI > yodals.txt
diff -u yodals.txt ${YODA_TESTS_SRC}/yodals-ref.txt

rm -f yodals.txt
