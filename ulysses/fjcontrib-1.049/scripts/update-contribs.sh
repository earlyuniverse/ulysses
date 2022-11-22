#!/bin/bash
#
# Usage:
#  update-contribs.sh
#
# update the list of all contribs.
# 
# This does the following
#
#  - run svn up to get the updates of all the scripts and latest contrib list
#
#  - build the list of contribs:
#     . by default, these are read from 'contribs.svn'
#     . local requests can be specified in 'contribs.local'
#       (in which case they take precedence)
#
#  - for all the contribs that have to be updated, check the currently
#    installed version and, if it does not match the one in
#    contribs.svn (or contribs.local), ask if the user wants to update
#    them. If yes, update the installed version.
#
# parse some command line options of the form
#   update-contribs.sh -<option>
#
# supported options at the moment are
#   -h      show a help message
#   --force assume the answer to every question is "yes"
if [[ $# -ge 1 && x"$1" == x'-h' ]]; then
    echo
    echo "Usage: "
    echo "       $0 [--force] [ContribName [version]] "
    echo 
    echo "- without any arguments, all contribs are updated (or downloaded if missing)"
    echo "- with the ContribName argument, just that contrib is updated"
    echo "- with additionally the version argument, the contrib is updated"
    echo "  (or switched) to the requested version. E.g. 'trunk' or 'tags/1.0' "
    echo
    exit 0
fi
default_yesno_answer=""
if [[ $# -ge 1 && x"$1" == x'--force' ]]; then
    echo "Assuming 'yes' as an answer to all questions"
    default_yesno_answer="yes"
fi
. `dirname $0`/internal/common.sh
    
internal_directories="_,scripts,Template,data,_"

#----------------------------------------------------------------------
# update svn

# get the revision of this file
script_current_version=$(svn info $0 | grep "^Revision: " | sed 's/Revision: //')

# perform the svn update
echo "-----------------------------"
echo "Updating top-level directory:"
echo "-----------------------------"
svn up || { echo "Failed to update svn. Aborting"; exit 1; }

# check if this script has been updated
script_new_version=$(svn info $0 | grep "^Last Changed Rev: " | sed 's/Last Changed Rev: //')
if [[ "$script_new_version" -gt "$script_current_version" ]]; then
    echo "update-contribs.sh has been updated. Re-running the new version."
    $0 || { exit 1;}
    exit 0
fi

#----------------------------------------------------------------------
# if there are two arguments, just call switch-to-version

if [[ $# -gt 1 ]]; then

    # just call switch-to-version
    `dirname $0`/internal/switch-to-version.sh $* || exit 1
    exit 0
fi

#----------------------------------------------------------------------
# update all the contribs in the contribs.svn file (plus the ones in contribs.local) or only the one
# specified through the command line
if [[ $# -gt 0 && x"$1" != x'--force' ]]; then
    contribs_list=$1
else
    contribs_list=$(cat contribs.svn | grep -v '^#' | grep -v '^$' | awk '{print " "$1" "}')
    # Check if a contribution mentioned in contribs.local was already present in contribs.svn
    # Only add the entry to contribs_list if it wasn't there before (the version number will
    # be dealt with later on)
    if [[ -e contribs.local ]]; then
      for contrib_in_local in `cat contribs.local | grep -v '^#' | grep -v '^$' | awk '{print $1}'`; do
	if [[ " $contribs_list " != *" $contrib_in_local "* ]]; then
	   contribs_list=$contribs_list" $contrib_in_local"
        fi
      done	
    fi		
fi 

echo
echo "--------------------------------------------"
echo "Checking for updates of individual contribs:"
echo "--------------------------------------------"
for contrib in $contribs_list; do
    # get the version numbers in contribs.svn file and also from the locally
    # checked out contributions
    get_contrib_version ${contrib} contribs.svn   version_svn
    get_contrib_version ${contrib} local_svn version_local
    get_contrib_version ${contrib} contribs.local version_mine

    echo
    echo -n "${contrib}: "
    
    # if thers is a line in contribs.local, use the version specified there
    # to supersede the one in contribs.svn (which could be implicitly [None],
    # i.e. a particular contribution could be not mentioned there)
    requested_tag="default"
    old_version=""
    if [[ "${version_mine}" != "[None]" ]]; then 
        old_version="  [Overriding $version_svn from contribs.svn]" 
        version_svn="$version_mine"
	requested_tag="requested"
    fi
    
    # check which situation we are in
    if [[ "${version_svn}" == "${version_local}" ]]; then
        # match: nothing to do
	if [[ "$version_svn" != "["*"]" ]]; then
	    echo -e "you already have the $requested_tag version (${version_svn}).\nRunning svn up"
	    cd $contrib
	    svn up
	    cd ..
	else 
	    echo "you already have the $requested_tag version (${version_svn})"
	fi	
    else
        #skip this particular contribution (flagged by at least a "-" in place of the version number)
        if [[ "${version_svn}" =~ ^-+ ]]; then echo "Skipped"; continue; fi
	    
	# mismatch: show the versions and decide what to do
	# according to the type of mismatch
	echo ""
	echo "    $requested_tag version: "${version_svn}$old_version
	echo "    installed version: "${version_local}
	if [[ "${version_local}" == "[None]" ]]; then
	    # the local version does not exist! Ask if we want to install it
	    #get_yesno_answer "  Do you want to install the $requested_tag version?" "$default_yesno_answer" || {
	    `dirname $0`/internal/switch-to-version.sh $contrib $version_svn || exit 1
	    #}
	elif [[ "${version_local}" == "[NoSVN]" ]]; then
	    echo "You have an unversionned copy of $contrib in the way. It will not be updated."
	else
	    # the local version exists! Ask if we want to update it
	    get_yesno_answer "  Switch from the installed version to the $requested_tag one?" "$default_yesno_answer" || {
		`dirname $0`/internal/switch-to-version.sh $contrib $version_svn || exit 1
	    }
	fi
        echo
    fi
done

#----------------------------------------------------------------------
# now do the opposite: for each local contrib, check if it exists in
# the supported lists
#
# Note that we discard any directory that does not point to a tagged
# version of a contrib

for contrib in $(ls -d */ || sed 's/\/*//g'); do
    # discard the fjcontrib dirs
    if [[ "$internal_directories" == *",${contrib},"* ]]; then
	continue
    fi

    get_svn_info $contrib mode version

    if [[ "$version" == "tags/"* ]]; then
	echo "${contrib}: your local copy ($version) does not appear in the default svn-suppoerted list."
	get_yesno_answer "  Do you want to remove the local version?" "$default_yesno_answer" || {
	    rm -Rf $contrib
	}
	echo
    fi
done
echo
