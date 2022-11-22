#!/bin/bash
#
# Usage:
#   scripts/new-contrib-from-template.sh <new_contrib_name>
#
# create the structure of a new contrib

#------------------------------------------------------------------------
# get the contrib name
if [ "x$1" == "x" ]; then
    echo ""
    echo "Usage:"
    echo "  $0 <new_contrib_name>"
    echo ""
    exit 1
fi

contrib=${1%/}

#------------------------------------------------------------------------
# make sure the name has not already been used
if [ -e $contrib ]; then
    echo "The $contrib contrib already exists. Please choose a different name"
    exit 1
fi

echo "Creating contrib "$1

contrib_lower=`echo ${contrib} | tr A-Z a-z`
contrib_upper=`echo ${contrib} | tr a-z A-Z`
date=`date "+%Y-%m-%d"`
user=`whoami`


#------------------------------------------------------------------------
# create the structure
mkdir $contrib

for fn in `dirname $0`/internal/Template/*; do
    stripped=${fn##*internal/}
    echo "  creating "${stripped//Template/${contrib}}
    sed "s/Template/${contrib}/g;\
         s/template/${contrib_lower}/g;\
         s/TEMPLATE/${contrib_upper}/g;\
         s/20XX-XX-XX/${date}/g;\
         s/xxxx@localhost/${user}@localhost/g"\
         ${fn} > ${stripped//Template/${contrib}}
done
echo "----------------------------------------------------------------------"
echo "$contrib successfully created from Template. Rerun ./configure"
echo "for it to be included in your builds."
echo
echo "Once you are ready to make it public, write to "
echo "fastjet@projects.hepforge.org "
echo "to obtain write access to the fastjet-contrib svn repository "
echo
echo "You may then start to upload your contrib by running "
echo
echo "    scripts/register-new-contrib.sh ${contrib}"
echo
echo "and following the instructions (details are in the README file)"
echo "----------------------------------------------------------------------"
