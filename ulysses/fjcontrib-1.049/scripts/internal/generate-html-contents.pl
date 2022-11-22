#!/usr/bin/perl -w
#
# Generate an html contents file for the current version of fjcontrib
#
# To be run from the main directory of a checkout of fjcontrib.
# Best run from a directory that is a clean checkout of a tag.
#
# set to 1 to sort contribs alphabetically, 0 otherwise
$sort=1;
# set to 1 to include release date taken from svn tags, 0 otherwise.
$include_date=1;

$versions="contribs.svn";

# get svn info from one common location
if (! exists $ENV{'svn_read'}) {
  print STDERR '$svn_read environment variable must be set, but is not (source common.sh first?)',"\n";
  exit(-1)
} else {
  $svn_read=$ENV{'svn_read'};
  print STDERR '$svn_read is: ', $svn_read,"\n";
}

$svn=$svn_read."/contribs/";
# link to browse the svn
$svnBrowse='https://phab.hepforge.org/source/fastjetsvn/browse/contrib/contribs/';
# this ensures that one doesn't get the whole "blame" info, which is ugly
$svnPost='?as=source&blame=off';

$topversion=`head -1 VERSION`;
chomp $topversion;

# read in contribs.svn file, fill contribs hash
open (VERSIONS, "<$versions") || die "Could not open $versions";
%contribs_hash=();
$contribs_array=();
while ($line = <VERSIONS>) {
  if ($line =~ /^\s*([a-z][^\s]*)\s+([^\s]*)/i) {
    $contrib = $1;
    $version = $2;
    push @contribs_array, $contrib;
    $contribs_hash{$contrib} = $version;
  }
}

# sort contribs by alphabetical order
if ($sort) { @contribs_array = sort keys %contribs_hash; }

# write out html table
$list='';
foreach ( @contribs_array ) { 
    $contrib = $_;
    $version = $contribs_hash{$_};
    if ($version =~ /^[0-9]/) {$version = "tags/$version";}
    ($textversion = $version) =~ s/tags\///;
    if($include_date) {
      # extract date of last version tag from svn
      print STDERR "getting date for $svn$contrib/$version\n";
      $date = `svn info $svn$contrib/$version | grep "Last Changed Date" | awk '{print \$4}'`;
      # One could  also use the --xml option and parse appropriately the output:
      # $date = `svn --xml list $svn$contrib/$version`;
      # At this stage, XML parsing is not implemented though.
      #
      #print $contrib." ".$date."\n";
    }
    $list .= "<tr> <td class=\"contribname\"> 
                   <a href=\"$svnBrowse$contrib/$version/\">$contrib</a>
               </td> <td style=\"{text-align:center;}\"> $textversion </td>";
    if ($include_date) {$list .= "<td>$date</td>";}
    $list .=  "<td>";
    if (-e "$contrib/README") {
      $list .= '<a href="'.$svnBrowse.$contrib.'/'.$version.'/README'.$svnPost.'">README</a> ';
    }
    if (-e "$contrib/NEWS") {
      $list .= '<a href="'.$svnBrowse.$contrib.'/'.$version.'/NEWS'.$svnPost.'">NEWS</a> ';
    }
    $list .= "</td></tr>\n";
}


$head='
<html>
<head>
<style type="text/css">
.contriblist table,.contriblist tr,.contriblist td,.contriblist th {
  text-align:center;
  border:1px solid white;
  padding:4px;
  padding-right:26px;
}
td.contribname {
  background-color:#eeeeee;
  text-align:left;
}
th.contriblist {
  text-align:center;
  padding:6px;
  padding-left:16px;
  padding-right:16px;
  background-color:#dddddd;
}
td.spanned {
  text-align:center;
  background-color:#dddddd;
}
</style>
</head>
<body>
Version '.$topversion.' of FastJet Contrib is distributed with the following packages<p>

<table class="contriblist">
<tr><th class="contriblist">Package</th> 
    <th class="contriblist">Version</th>';
if($include_date) {$head .= '
    <th class="contriblist">Release date</th>';}
$head .= '
    <th class="contriblist">Information</th> </tr> 
';

$tail='
</table>


</body>
</html>
';

print $head, $list, $tail;
