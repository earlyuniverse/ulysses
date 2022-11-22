#! /usr/bin/env bash

yoda2root ${YODA_TESTS_SRC}/test1.yoda
ST=$?
if [[ $ST = 2 ]]; then exit 0; else exit $ST; fi #< Failure to find the ROOT module is "ok"
