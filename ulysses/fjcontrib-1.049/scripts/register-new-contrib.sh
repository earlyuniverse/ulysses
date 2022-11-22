#!/bin/bash
#
# Usage:
#  register-new-contrib.sh <ContribName>
#
# uploads a new contrib to the svn. Note that only the svn directory
# structure will be created and the files in the contrib's directory
# will be simply moved locally but not added to the svn repository.
#
# At the end of the procedure, the Contrib's directory will point to
# the trunk

# load common setup
. `dirname $0`/internal/common.sh

# make sure we have an argument
contrib=${1%/}
if [ -z $contrib ]; then
    echo "Usage:"
    echo "  register-new-contrib.sh <ContribName>"
    echo "A contrib name has to be specified"
    exit 1
fi

echo
echo "This script will create the directory structures for $contrib in the svn repository"
echo "It will also arrange for the $contrib directory to have a .svn/ pointing to "
echo
echo  "  " $svn_write/contribs/$contrib/trunk
echo
echo "so that you can perform svn operations on your contrib"
echo
get_yesno_answer "Do you wish to continue?" && exit 1
echo

# check that the contrib name exist locally
echo "  performing sanity checks"
if [ ! -d $contrib ]; then
    echo "  $contrib does not exist locally"
    exit 1
fi


# check that the current contrib does not exist
if [ ! -z "`svn ls $svn_read/contribs | grep '^'$contrib'/$'`" ]; then
    echo "  $contrib is the name of an already-existing contrib. Aborting..."
    exit 1
fi

#TODO check that the directory contains "expected" files

# create the svn structure for that contrib
echo "  Creating the svn folder structure (a password may be requested)"
svn mkdir -m "Creating the basic svn structure for contrib $contrib" \
    $svn_write/contribs/$contrib \
    $svn_write/contribs/$contrib/trunk \
    $svn_write/contribs/$contrib/tags \
    $svn_write/contribs/$contrib/branches \
 || { echo "Failed to create the svn directory structure. Aborting."; exit 1;}

# now we need to have "contrib" point to that svn location.
#  - move the existing one out of the way
echo "  Moving $contrib to ${contrib}.bak"
mv $contrib ${contrib}.bak

#  - make a checkout
echo "  Checking out the svn trunk for $config (a password may be requested)"
svn co ${svn_write}/contribs/$contrib/trunk $contrib

# and finally, add all the contrib files to the svn (with some trickery to make
# sure dot files are moved by mv * target/
default_dotglob_status="-u"
if [ -z "`shopt dotglob | grep "off"`" ]; then
    default_dotglob_status="-s"
fi
shopt -s dotglob
echo "  Copying all files from ${contrib}.bak into the new $contrib/ directory"
cp -a ${contrib}.bak/* ${contrib}
shopt $default_dotglob_status dotglob

echo
echo "You should now have a $contrib/ directory, registered under svn and with your original files"
echo
echo "Don't forget to enter the $contrib/ directory, and run svn add, svn commit"
echo "to populate the repository"
