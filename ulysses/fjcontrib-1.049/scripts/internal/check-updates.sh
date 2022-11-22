#!/bin/bash
#
# check if the versions in contrib.svn correspond to the latest tags

# preamble
if [ x`which tput` != "x" ]; then
   GREEN=$(tput setaf 2)
   RED=$(tput setaf 1)
   NORMAL=$(tput sgr0)
fi
. `dirname $0`/common.sh

# get the list of contribs (discard "graveyard"
contrib_list=`svn ls $svn_read/contribs/ | grep -v graveyard | sed 's/\///g'`

# loop over contribs
printf "  %-25s %-15s %-15s\n" "contrib" "contribs.svn" "svn tag"
printf "  %-25s %-15s %-15s\n" "-------" "------------" "-------"
for contrib in $contrib_list; do
    # check version in contribs.svn
    get_contrib_version $contrib contribs.svn version_included
    version_included=`echo $version_included | sed 's/.*\///'`
    
    # check latest svn tag
    version_tag=`svn ls $svn_read/contribs/${contrib}/tags | grep -E "^[0-9].[0-9].[0-9]/$" | tail -n1 | sed 's/\///g'`

    # see if that agrees
    if [ x"$version_included" == x"$version_tag" ]; then
        col=$GREEN
    elif [ x"$version_included" == x"[None]" ]; then
        if [ x"$version_tag" == x"" ]; then
            col=$NORMAL
        else
            col=$RED
        fi
    elif [ x"$version_included" \> x"$version_tag" ]; then
        col=$NORMAL
    else
        col=$RED
    fi
    printf "%s  %-25s %-15s %-15s%s\n" ${col} $contrib $version_included $version_tag $NORMAL
done




