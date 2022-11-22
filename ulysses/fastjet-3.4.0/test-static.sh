#!/bin/bash
# Script for testing whether everything works with the old (pre 2.3) build 
# framework. This compiles things and runs a few tests.

# Usage
# ./test-script.sh [-c] [-j2]

# following was designed for when pxcone was part of test script.
# # set G77LIBDIR if the user has not set it in its environment
# if [[ ! $G77LIBDIR ]] ; then
#     dir=`locate libg2c.a | head -1`
#     dir=`dirname $dir` 
#     echo Setting G77LIBDIR to $dir
#     echo Press return to accept, or input a different directory
#     read newdir
#     if [[ $newdir == "" ]]; then
#         export G77LIBDIR=$dir
#     else
#         export G77LIBDIR=$newdir
#     fi
# fi

# set CLEAN to "empty string" in order NOT to clean
CLEAN=""
MAKEARGS="-f makefile.static "

# flag -c in input to clean, -jN to make
for arg in $*
do
  if [[ $arg == "-c" ]]; then
      CLEAN="1"
  elif [[ $arg == "-j2" ]]; then
      MAKEARGS=$MAKEARGS" "$arg
  elif [[ $arg == "-j4" ]]; then
      MAKEARGS=$MAKEARGS" "$arg
  elif [[ $arg == "-j8" ]]; then
      MAKEARGS=$MAKEARGS" "$arg
  fi
done 


# output file
OUTPUT=/tmp/output$$

# build the main package
pushd src
if [ $CLEAN ]; then make $MAKEARGS clean || exit -1; fi
make $MAKEARGS install
popd

# build the plugins
pushd plugins
if [ $CLEAN ]; then make $MAKEARGS clean || exit -1; fi
make $MAKEARGS || exit -1
popd

# build the examples -- tests will be postponed to later
pushd example
if [ $CLEAN ]; then make $MAKEARGS clean || exit -1; fi
make $MAKEARGS fastjet_example  || exit -1
## (./fastjet_example < data/single-event.dat | tee $OUTPUT) || exit -1
make $MAKEARGS fastjet_areas || exit -1
## (./fastjet_areas < data/single-event.dat | tee -a $OUTPUT) || exit -1
make $MAKEARGS fastjet_example_v1_interface || exit -1
make $MAKEARGS fastjet_subtraction || exit -1
make $MAKEARGS fastjet_timing || exit -1
#make $MAKEARGS fastjet_timing_plugins || exit -1
popd

# build some plugin examples
pushd plugins/usage_examples
if [ $CLEAN ]; then make $MAKEARGS clean || exit -1; fi
#make $MAKEARGS cdfmidpoint_example many_algs_example pxcone_example siscone_example || exit -1
make $MAKEARGS cdfmidpoint_example many_algs_example siscone_example || exit -1
## (./many_algs_example < data/single-event.dat | tee -a $OUTPUT) || exit -1
popd

# now run tests from the tests-compare.sh
export srcdir=.
./test-compare.sh


## # mv the output file to a "decent" place (remove 
## #mv $OUTPUT  test-script-output.txt
## grep -v 'CGAL' $OUTPUT > test-script-output.txt
## echo 
## echo -----------------------------------------------------------
## echo "Comparing output from these runs (test-script-output.txt) "
## echo "to the expected output (test-script-output-orig.txt)"
## echo -----------------------------------------------------------
## grep -v -e '#' -e 'SISCone' test-script-output.txt > output1.tmp
## grep -v -e '#' -e 'SISCone' test-script-output-orig.txt > output2.tmp
## #diff --ignore-matching-lines="#" --ignore-matching-lines=SISCone test-script-output.txt  test-script-output-orig.txt > $OUTPUT
## #DIFF=`diff --ignore-matching-lines="#" --ignore-matching-lines=SISCone test-script-output.txt  test-script-output-orig.txt`
## diff output1.tmp output2.tmp > $OUTPUT
## DIFF=`diff output1.tmp output2.tmp`
## rm output1.tmp output2.tmp
## if [[ $DIFF ]]; then cat $OUTPUT
## else 
## echo Results are identical
## fi
## 
## # remove temporary output file
## rm -f $OUTPUT
