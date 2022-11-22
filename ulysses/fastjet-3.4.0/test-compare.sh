#!/bin/bash

# This file aims at comparing the output of the different algorithms 
# with the benchmark results stored in test-algorithms-orig.dat
#
# To do that, we'll follow the following recipe:
#   - get the list of available algorithms
#     Note that at the same time, we'll build the info to clear the 
#     original output later on
#   - run the various algorithms using fastjet_timing_plugins
#     send the output to a tmp file
#   - clear orig and tmp files
#   - compare them
# 
# Specific notes:
#   - Tested algorithms: 
#       kt, cam, antikt, ee_kt, ee_cam, ee_antikt,
#       siscone, siscone_spheri, jetclu, midpoint, pxcone,
#       d0runiicone, trackjet
#   - What needs to be cleared in the orig file:
#       . the lines starting with "<alg>:" for the algorithms
#         that won't be used
#       . the lines starting with #
#       . the lines containing "version" (version number)
#       . the lines containing "strategy" (strategy night differ)
#       . the lines containing ":#" (comments for the used algs)
#       . the lines containing "SISCone" (version number)
#       . the lines containing "CGAL" (not necessarily available)
#   - What needs to be cleared in the test file:
#       . all comment lined starting with '#'
#       . the lines containing "version" (version number)
#       . the lines containing "strategy" (strategy night differ)
#       . the lines containing "SISCone" (version number)
#       . the lines containing "CGAL" (not necessarily available)
#   - For native algs, we'll fix the strategy so that the output 
#     matches (Best could be dangerous in case CGAL gets used)
#   - additional things tested: the output of fastjet_example and
#     fastjet_areas that will be placed at the beginning of the file
#   - for algs that are not in the plugin <Alg>Plugin, specify an 
#     additional tag using
#        <alg>:<plugin_name>
#   - additional parameters for 'fastjet_timing_plugins' can be 
#     specified using ',' to separate them
#   - another event file can be specified using @<eventfile>
#     important: only the name of teh file has to be specified
#                and it has to be in example/data/
# 


# if srcdir is not defined, set it to .
# this allows to run outside of make check
if test -z ${srcdir}; then
    echo "setting srcdir to ."
    srcdir="."
fi



#======================================================================
# first check that the windows include file has the correct
# version number. Strictly speaking not something that 
# should be done in make check, but should help guarantee
# that we do not release a copy with the wrong windows version 
# number
packname=`grep '^ *AC_INIT' ${srcdir}/configure.ac | sed -e 's/AC_INIT(//' -e 's/\[//g' -e 's/\]//g' -e 's/)//'`
packver=`echo $packname | sed 's/.*,//g'`
winver=`grep FASTJET_PACKAGE_VERSION ${srcdir}/include/fastjet/config_win.h | sed 's/.*VERSION *"//' | sed 's/"//'`
if [[ $winver != $packver ]]; then
  echo "ERROR: config_win.h version number not compatible with true version number"
  exit 1
fi

#======================================================================
# check that 
#    fastjet-config --cxxflags --libs
# returns the same as
#    fastjet-config --cxxflags; fastjet-config --libs
arguments_full=`./fastjet-config --cxxflags --libs | sed 's/\n/ /'`
arguments_flags=`./fastjet-config --cxxflags | sed 's/\n/ /'`
arguments_libs=`./fastjet-config --libs | sed 's/\n/ /'`
arguments_recon=${arguments_flags}" "${arguments_libs}
if [[ "$arguments_full" != "$arguments_recon" ]]; then
  echo "ERROR: fastjet-config --cxxflags --libs does not combine fastjet-config --cxxflags and fastjet-config --libs"
  echo "  fastjet-config --cxxflags       : "$arguments_flags
  echo "  fastjet-config --libs           : "$arguments_libs
  echo "  fastjet-config --cxxflags --libs: "$arguments_full
  exit 1
fi


#======================================================================
# now run the real tests
# first build the list of algs to run
echo -----------------------------------------------------------
echo "Checking which algorithms are available for testing"
echo -----------------------------------------------------------
tested_algs="kt cam antikt genkt,1.0 genkt,0.0 genkt,-1.0 eekt,-excld,2.0@single-ee-event.dat eegenkt,1.0@single-ee-event.dat eegenkt,0.0@single-ee-event.dat eegenkt,-1.0@single-ee-event.dat"
untested_algs=""

Rvalues="0.4 0.7 1.0"
extra_args="-incl 0.0"

echo ":#"          >  clear_patterns.orig
echo "^#"          >> clear_patterns.orig
echo "version"     >> clear_patterns.orig
echo "[Ss]trategy" >> clear_patterns.orig
echo "CGAL"        >> clear_patterns.orig
echo "SISCone"     >> clear_patterns.orig  # avoids problems w version numbers
echo "pxcone: +[a-zA-Z*]" >> clear_patterns.orig   # special treatment for PxCone whose fortran output
echo "pxcone: *$"         >> clear_patterns.orig   # occurs in non-predicatble position (flushing issue)
echo "WARNING"            >> clear_patterns.orig   # occurs in non-predicatble position (flushing issue)
cp clear_patterns.orig clear_patterns.tmp

# note: algs specified as alg:name mean that 'name' has to be checked for the 
#       availability of 'alg'
for plugin_tag in siscone sisconespheri:siscone jetclu:cdfcones midpoint:cdfcones pxcone d0runiicone trackjet atlascone cmsiterativecone eecambridge@single-ee-event.dat:eecambridge eecambridge,-ycut,0.01@single-ee-event.dat:eecambridge jade,-excly,0.01@single-ee-event.dat:jade d0runipre96cone:d0runicone d0runicone:d0runicone gridjet ; do
    plugin=${plugin_tag%%:*}
    tag=${plugin_tag##*:}

    tag_upper=`echo ${tag} | tr a-z A-Z`

    if [[ -n `grep "define FASTJET_ENABLE_PLUGIN_${tag_upper}" include/fastjet/config_auto.h` ]]; then
	tested_algs=${tested_algs}" "${plugin}
    else
	untested_algs=${untested_algs}" "${plugin}
	echo "^${plugin}:" >> clear_patterns.orig
    fi
done

# build the output to be compared wityh the original one
## for regenerating the orig output: echo "blahblahthiswillneverhappen" > clear_patterns.tmp
echo -----------------------------------------------------------
echo "Running 'fastjet_example < data/single_event.dat'"
echo -----------------------------------------------------------
example/fastjet_example < ${srcdir}/example/data/single-event.dat | grep -v -E -f clear_patterns.tmp > output.tmp

echo
echo -----------------------------------------------------------
echo "Running 'fastjet_areas < data/single_event.dat'"
echo -----------------------------------------------------------
example/fastjet_areas < ${srcdir}/example/data/single-event.dat | grep -v -E -f clear_patterns.tmp >> output.tmp

# run the algorithms to be tested
echo -----------------------------------------------------------
echo "Running 'fastjet_timing_plugins "${extra_args}" < data/single_event.dat' on all algs"
echo "  tested  : "${tested_algs}
echo "  untested:" ${untested_algs}
echo "  R values: "${Rvalues}
echo -----------------------------------------------------------
for alg in ${tested_algs}; do
    # check if we have to run on a separate event
    # Note that this is specified using @<event_file>
    alg_part=${alg%%@*}
    if [ "${#alg}" -eq "${#alg_part}" ]; then
	event_part="single-event.dat"
    else
	event_part=${alg##*@}
    fi

    # additional parameters for 'fastjet_timing_plugins' can be specified using ',' to separate them.
    # we thus have to replace ',' by ' ' when we run fastjet_timing_plugins
    for R in ${Rvalues}; do
	example/fastjet_timing_plugins -${alg_part//,/ } ${extra_args} -r ${R} < ${srcdir}/example/data/${event_part} 2>&1 \
	    | awk "{print \"${alg}:\"\$0}"  \
	    | grep -v -E -f clear_patterns.tmp \
            >> output.tmp 
    ##      | awk "{if (\$2 == \"exclusive\"){ exit;}; print \"${alg}:\"\$0}"  >> output.tmp
    done
done
## for regenerating the orig output: cp output.tmp test-script-output-orig.txt

echo
echo -----------------------------------------------------------
echo "Comparing output from these runs (test-script-output.txt) "
echo "to the expected output (test-script-output-orig.txt)"
echo -----------------------------------------------------------
# clear the original output for comment lines and untested algorithms
grep -v -E -f clear_patterns.orig  ${srcdir}/test-script-output-orig.txt >  output_orig.tmp

# 4. perform the diff
DIFF=`diff output.tmp output_orig.tmp`
diff output.tmp output_orig.tmp > test-script-output.tmp

# 5. show result
echo "Tested plugins: "${tested_plugins}
if [[ -n $DIFF ]]; then 
  cat test-script-output.tmp
  exit 1;
else
  echo Results are identical

  rm test-script-output.tmp
  rm clear_patterns.orig
  rm clear_patterns.tmp
  rm output.tmp output_orig.tmp
fi
