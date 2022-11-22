
#!/bin/bash

# This script, meant for "internal use", packages all the
# latest version of the contributions into a tar.gz file,
# either reading the list from contribs.svn, or
# exploiting the list returned by running "./configure --list"
# in the top directory (see commented line)
#
#  Usage: run ./scripts/internal/package.sh from the top directory
#

version=`cat VERSION`
packagename=fjcontrib-$version
files="AUTHORS COPYING Makefile.in README VERSION NEWS configure data ChangeLog utils/check.sh utils/install-sh"
# get contribs list from contribs.svn
files="$files "$(cat contribs.svn | grep -v "^#" | grep -v "^$" | awk '{printf $1" "}')
# get contribs list from configure (may include "unofficial" local contribs)
#files="$files "`./configure --list`

echo
echo Making a distclean, so as to only package relevant files
echo
make distclean
echo
echo Creating version $version
echo Including: $files
tar --exclude=".svn"  --transform "s,^,$packagename/," -zcf $packagename.tar.gz $files

exit


