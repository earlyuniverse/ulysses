#!/bin/bash

set -e

yodascale -c '.* 90.0' ${YODA_TESTS_SRC}/test2.yoda
mv test2-scaled.yoda scale1.yoda

yodascale -c '.* 10x'  ${YODA_TESTS_SRC}/test2.yoda
mv test2-scaled.yoda scale2.yoda

yodadiff scale1.yoda scale2.yoda

rm -f scale1.yoda scale2.yoda

# yodascale -c '.* 0.25x' ${YODA_TESTS_SRC}/test3.yoda
# mv test3-scaled.yoda scale1.yoda

# yodascale -c '.* 111.0' ${YODA_TESTS_SRC}/test3.yoda
# mv test3-scaled.yoda scale2.yoda

# yodadiff scale1.yoda scale2.yoda

# rm -f scale1.yoda scale2.yoda
