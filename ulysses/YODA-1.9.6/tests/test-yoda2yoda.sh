#!/bin/bash

set -e

yoda2yoda ${YODA_TESTS_SRC}/rivetexample.yoda y2y_1.yoda -m DELPHI
yodacnv   ${YODA_TESTS_SRC}/rivetexample.yoda y2y_2.yoda -m DELPHI

yodadiff y2y_1.yoda y2y_2.yoda

yoda2yoda --as-scatters ${YODA_TESTS_SRC}/rivetexample.yoda y2y_3.yoda -m DELPHI

# Check that every analysis object is a scatter
all_aos=$(yodals y2y_3.yoda | tail -n +2 | awk '{print $2}' | wc -l)
scat_aos=$(yodals y2y_3.yoda | tail -n +2 | awk '{print $2}' | grep 'Scatter' | wc -l)
[[ ${all_aos} -eq ${scat_aos} ]]

rm -f y2y_1.yoda y2y_2.yoda
