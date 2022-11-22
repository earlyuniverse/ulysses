#!/bin/bash
#
# make a full release of the current trunk
# Then produce a tarball

# include function and svn location definitions, etc.
. `dirname $0`/common.sh

dry_run=0
if [[ "$1" == "--dry-run" ]]; then
    dry_run=1
    echo "--------------------------------------------------"
    echo "                   DRY RUN                        "
    echo "--------------------------------------------------"
fi

#========================================================================
# svn sanity checks
#========================================================================

# make sure that everything is committed
echo
echo "Checking for pending modifications or updates (this may take a few seconds...)"
if [[ ! -z  "`svn status --show-updates | grep -v "^?" | grep -v "^Status"`" ]]; then
    echo
    echo "WARNING: There are pending modifications or updates:"
    echo
    svn status --show-updates | grep -v '^\?'
    echo
    get_yesno_answer "Are you really sure you want to proceed?" &&  exit 1
    echo
else
    echo "All files are up to date relative to the svn"
fi

# make sure there is a VERSION and it does not already exist
version=`head -n1 VERSION`
if [[ ! -z $(svn ls $svn_read/tags | grep "^$version/") ]]; then
    echo "Version $version of fjcontrib already exists. Aborting"
    exit 1
fi

echo
echo "The contribs.svn file points to the following contrib versions"
echo
echo "----------------------------------------------------------------"
grep -v '^#' contribs.svn
echo "----------------------------------------------------------------"
echo

get_yesno_answer "Do you want to proceed with the release of fjcontrib-$version?" &&  exit 1


#========================================================================
# check that the tools in contribs.svn behave OK
#========================================================================
# get a clean checkout to perform sanity checks
svn co $svn_read/trunk fjcontrib-$version || { echo "Failed to do the svn checkout"; exit 1; }
cd fjcontrib-$version
echo "------------------------------------------------------------------------"
echo "Getting the contribs"
echo "------------------------------------------------------------------------"
if ./scripts/update-contribs.sh --force; then
    echo "Success."
    echo
else
    echo "Failed."
    echo
    cd ..
    exit 1
fi

echo "------------------------------------------------------------------------"
echo "Configuring"
echo "------------------------------------------------------------------------"

# we need to determine whether to use fastjet-config from the path or
# use the one from the configure invocation in the trunk
is_in_path="yes"
which fastjet-config > /dev/null || is_in_path="no"

trunk_version=""
if [[ -e "../Makefile" ]]; then
    trunk_version=$(head -n3 ../Makefile | tail -n1 | grep "\--fastjet-config=" | sed 's/.*--fastjet-config=//;s/ .*$//')
fi

if [[ -z "$trunk_version" ]]; then
    if [[ "$is_in_path" == "no" ]]; then
	echo "fastjet-config is not in your path and cannot be obtained from the trunk configuration. Aborting."
	cd ..
	exit 1
    else
	echo "Using fastjet-config from your path"
	configure_options=""
    fi
else
    if [[ "$is_in_path" == "no" ]]; then
	echo "using fastjet-config from the trunk configuration"
	configure_options=" --fastjet-config=${trunk_version}"
    else
	echo "fastjet-config can be either taken from your path or from $trunk_version."
	configure_options=""
	get_yesno_answer "Do you want to use the one from your trunk?" || {
	    configure_options=" --fastjet-config=${trunk_version}"
	}
	    
    fi
fi

if ./configure $configure_options; then
    echo "Success."
    echo
else
    echo "Failed."
    echo
    cd ..
    exit 1
fi

echo "------------------------------------------------------------------------"
echo "Running make check"
echo "------------------------------------------------------------------------"
if make -j4 check; then
    echo "Success."
    echo
else
    echo "Failed."
    echo
    cd ..
    exit 1
fi

echo "------------------------------------------------------------------------"
echo "Running make fragile-shared"
echo "------------------------------------------------------------------------"
if make -j4 fragile-shared; then
    echo "Success."
    echo
else
    echo "Failed."
    echo
    cd ..
    exit 1
fi

cd ..
rm -Rf fjcontrib-$version
if [ -d fjcontrib-$version ]; then
    echo "fjcontrib-$version still present. Aborting"
fi

#========================================================================
# tag the release
#=======================================================================
if (( ${dry_run} )); then
    echo "Dry run: skipping the release tag"
else
    echo
    get_yesno_answer "Confirm you want to tag the release and make a tarball?" &&  exit 1
    echo
    
    echo "------------------------------------------------------------------------"
    echo "Making a tag of fjcontrib version $version"
    echo "------------------------------------------------------------------------"
    echo svn copy -m "tagging fjcontrib-$version" $svn_write/trunk $svn_write/tags/$version
         svn copy -m "tagging fjcontrib-$version" $svn_write/trunk $svn_write/tags/$version
fi

#========================================================================
# produce a tarball
#========================================================================
if (( ${dry_run} )); then
    echo "------------------------------------------------------------------------"
    echo "Dry run: checking out the trunk build the fjcontrib tarball"
    echo "------------------------------------------------------------------------"
    # using svn_write, because the http access sometimes doesn't
    # immediately see the up to date svn repository(?!)
    echo svn co $svn_read/trunk fjcontrib-$version
    svn co $svn_read/trunk fjcontrib-$version || { echo "Failed to checkout the fjcontrib trunk"; exit 1; }
else
    echo "------------------------------------------------------------------------"
    echo "Checking out tags/$version of fjcontrib"
    echo "------------------------------------------------------------------------"
    # using svn_write, because the http access sometimes doesn't
    # immediately see the up to date svn repository(?!)
    echo svn co $svn_write/tags/$version fjcontrib-$version
    svn co $svn_write/tags/$version fjcontrib-$version || { echo "Failed to checkout the new released version tags/$version"; exit 1; }
fi
cd fjcontrib-$version
echo

echo "------------------------------------------------------------------------"
echo "Getting the contribs"
echo "------------------------------------------------------------------------"
if ./scripts/update-contribs.sh --force; then
    echo "Success."
    echo
else
    echo "Failed."
    echo
    cd ..
    exit 1
fi

# # get rid of a few things for developers and "svn-users" only
# mkdir tmp
# for fn in check.sh install-sh; do
#     mv scripts/internal/${fn} ./tmp
# done
# rm -Rf scripts
# mkdir scripts
# mkdir scripts/internal
# for fn in tmp/*; do
#     mv $fn scripts/internal/${fn#tmp/}
# done
# rm DEVEL-GUIDELINES

cd ..
echo "------------------------------------------------------------------------"
echo "Producing fjcontrib-$version.tar.gz"
echo "------------------------------------------------------------------------"
tar --exclude=".svn" \
    --exclude="fjcontrib-$version/contribs.svn" \
  -czf fjcontrib-$version.tar.gz fjcontrib-$version
rm -Rf fjcontrib-$version
echo
echo "Success."
echo

#========================================================================
# update things on HepForge
#========================================================================
if (( ${dry_run} )); then
    echo "Dry run: not updating HEPForge"
    echo
    echo "Done"
    echo
    exit 0
fi

echo
get_yesno_answer "Confirm you want to upload to hepforge?" &&  exit 1
echo
echo "------------------------------------------------------------------------"
echo "Uploading to HepForge"
echo "------------------------------------------------------------------------"

echo "Uploading fjcontrib-$version.tar.gz"
scp fjcontrib-$version.tar.gz login.hepforge.org:$fastjet_web_dir/contrib/downloads/

mkdir hepforge_tmp
echo "Generating info for the webpage"
echo -n "$version" > hepforge_tmp/fjcversion.php
`dirname $0`/generate-html-contents.pl > hepforge_tmp/contents-$version.html
reldate=`date +"%e %B %Y"`
echo -n $reldate  > hepforge_tmp/fjcreldate.php

echo "Uploading info for the webpage"
scp hepforge_tmp/fjcversion.php hepforge_tmp/fjcreldate.php login.hepforge.org:$fastjet_web_dir/contrib/
scp hepforge_tmp/contents-$version.html login.hepforge.org:$fastjet_web_dir/contrib/contents/$version.html


echo "Ensuring fastjet group write access for new files on hepforge"
# the following is needed because group sticky bit is erroneously not set
# on the fastjet downloads directory, so group does not get set to fastjet
#ssh login.hepforge.org chgrp fastjet "~fastjet/downloads/fjcontrib-$version.tar.gz"
# now give fastjet group write permission on these files
ssh login.hepforge.org chmod g+w "$fastjet_web_dir/contrib/fjcversion.php" "$fastjet_web_dir/contrib/fjcreldate.php" "$fastjet_web_dir/contrib/contents/$version.html" "$fastjet_web_dir/contrib/downloads/fjcontrib-$version.tar.gz"
rm -Rf hepforge_tmp
echo
echo "Done"
echo


