#!/bin/bash

set -e

yodamerge ${YODA_TESTS_SRC}/test1.yoda ${YODA_TESTS_SRC}/test2.yoda -o merged12.yoda
yodadiff merged12.yoda ${YODA_TESTS_SRC}/merged12-ref.yoda

yodastack ${YODA_TESTS_SRC}/test1.yoda:2 ${YODA_TESTS_SRC}/test2.yoda:3.142 -o merged12pi.yoda
yodadiff merged12pi.yoda ${YODA_TESTS_SRC}/merged12pi-ref.yoda

rm -f merged12.yoda merged12pi.yoda
