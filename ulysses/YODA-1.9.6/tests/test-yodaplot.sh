#!/bin/bash

set -e

yodaplot -f svg -E MPL ${YODA_TESTS_SRC}/test1.yoda

rm -f foo_bar_baz.svg
