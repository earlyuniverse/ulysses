#!/bin/bash
#
# Usage:
#   switch_to_version ContribName Version
#
# switch the given contrib to the given version.
#
# In both cases, if the version is "trunk", we give it a svn_write
# access; otherwise, we use a read-only (http) access.

# load some common utilities
. `dirname $0`/common.sh

# make sure a contrib is specified
contrib=${1%/}
if [[ -z "$contrib" ]]; then
    echo "A contrib name has to be specified"
    exit 1
fi

# decide wich version number to use
version=$2
if [[ -z "$version" ]]; then
    # get the version from the svn list
    get_contrib_version $contrib contribs.svn version_svn
    if [[ "$version_svn" == "[None]" ]]; then
	echo "$contrib is not listed in the svn-supported list of contribs (contribs.svn)."
	echo "You have to specify the version yourself."
    	if [ ! -z "`svn ls $svn_read/contribs | grep '^'$contrib'/$'`" ]; then
    		echo "Possible tagged (released) versions are:"
    		svn ls $svn_read/contribs/$contrib/tags | sed 's/\/$//g;s/^/  tags\//g'
    	fi
	exit 1
    fi
    echo "  Using version $version_svn"
    version=$version_svn
fi

# decide which access mode we want
mode="ro"
if [[ "$version" == "trunk" ]]; then
    mode="rw"
fi
if [[ "${version}" =~ ^branches ]]; then
    mode="rw"
fi    

# if the version starts with a number, prefix it by "tags/"
if [[ "$version" =~ ^[0-9] ]]; then
    version="tags/$version"
fi

# get the current situation
get_svn_info $contrib current_mode current_version
if [[ "$version" == "$current_version" ]]; then
    if [[ "$mode" == "$current_mode" ]]; then
	if [[ "$version" != "["*"]" ]]; then
	    echo "  Already at the requested version. Simply running svn up"
	    cd $contrib
	    svn up
	    cd ..
	else 
	    echo "  Already at the requested version."
	fi	
	exit 0
    fi
fi

# make sure that the version exists
if [ -z "`svn ls $svn_read/contribs | grep '^'$contrib'/$'`" ]; then
    echo "${contrib} does not appear to be a valid contrib in the svn repository"
    exit 1
fi
if [[ "${version}" != "trunk" ]]; then
    # check if we're requesting a tag or a branch
    if [[ "$version" =~ ^branches ]]; then
        # we deal with a branch
        branch_version=${version#*/}
        if [ -z "`svn ls $svn_read/contribs/$contrib/branches | grep '^'$branch_version'/$'`" ]; then
	    echo "Version $version of $contrib does not exist. Exiting."
	    exit 1
        fi
    else
        # we deal with a tag
        tagged_version=${version#*/}
        if [ -z "`svn ls $svn_read/contribs/$contrib/tags | grep '^'$tagged_version'/$'`" ]; then
	    echo "Version $version of $contrib does not exist. Exiting."
	    exit 1
        fi
    fi
fi


# if the directory does not yet exist, check it out
if [[ "$current_version" == "[None]" ]]; then
    echo "  Checking out version ${version} of ${contrib}"
    if [[ "${mode}" == "rw" ]]; then
        svn co ${svn_write}/contribs/${contrib}/${version} $contrib
    else
        svn co ${svn_read}/contribs/${contrib}/${version} $contrib
    fi
    exit 0
fi

# if the directory exists but is not on svn, print an error
# TODO: allow to overwrite your local copy
if [[ "$current_version" == "[NoSVN]" ]]; then
    echo "You appear to have an unversionned copy of ${contrib}."
    echo "Please move it out of the way before updating ${contrib} to a versioned version."
    exit 1
fi

# check that there is no pending 
check_pending_modifications $contrib || {
    get_yesno_answer "Your local copy has modifications. Do you want to proceed with the update?" && {
	echo "Aborting."
	exit 1
    }
}

# move to the contrib's directory
cd $contrib

# in principle, we should use an "svn switch" here (change
# a directory in a given repos). But that does not allow
# to change the access type. So if the access type
# changes, we combine that with a svn relocate
#
# Options include:
#  - a version upgrade: svn switch
#  - trunk->stable: svn relocate + svn switch
#  - stable->trunk: svn switch + svn relocate
# Se we need to check if we are currently using a svn+ssh
# or http access

# if we're coming from the trunk, first switch to a read-only access
if [[ "$current_mode" == "rw" ]]; then
    echo "  (first changing the access-mode to read-only)"
    svn switch --relocate ${svn_write}/contribs/$contrib/${current_version} ${svn_read}/contribs/$contrib/${current_version}
    svn up
fi

# now switch to the correct location
svn switch ${svn_read}/contribs/$contrib/$version
# grant write access when needed
if [[ "${mode}" == "rw" ]]; then
    svn switch --relocate ${svn_read}/contribs/$contrib/${version} ${svn_write}/contribs/$contrib/${version}
    svn up
fi

cd ..
