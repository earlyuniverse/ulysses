# a list of definitions and tools that we'd like to ave an easy access
# to

# svn repositories for read and write access
#svn_read=http://fastjet.hepforge.org/svn/contrib
#svn_write=https://fastjet.hepforge.org/svn/contrib
#svn_write=svn+ssh://svn.hepforge.org/hepforge/svn/fastjet/contrib
#fastjet_web_dir=~fastjet/public_html
#svn_read=svn+ssh://vcs@phab.hepforge.org/source/fastjetsvn/contrib
svn_read=https://svn.hepforge.org/fastjetsvn/contrib
svn_write=svn+ssh://vcs@phab.hepforge.org/source/fastjetsvn/contrib
fastjet_web_dir=/hepforge/projects/fastjet/public_html

export svn_read svn_write fastjet_web_dir

# get the svn URL and fill 
#  - mode : ro if http:// access; rw otherwise
#  - version : the version info
#               [None]  is returned if the directory does not exist
#               [NoSVN] is returned if the directory is not under svn
#
#   get_svn_info  contrib  mode  version
function get_svn_info(){
    local __modevar=$2
    local __versionvar=$3

    # check if the directory exists
    if [[ ! -d $1 ]]; then
	eval $__modevar="[None]"
	eval $__versionvar="[None]"
	return 0
    fi

    cd $1

    # check if this is in svn
    svn info > /dev/null 2>&1 || {
	eval $__modevar="[NoSVN]"
	eval $__versionvar="[NoSVN]"
	cd ..
	return 0
    }
	
    # get the full URL
    svn_url=$(svn info | grep "^URL:" | sed 's/^URL: //')
    eval $__versionvar="${svn_url#*/$1/}"

    if [[ "$svn_url" == "http:"* ]]; then
	eval $__modevar="ro"
    else
	eval $__modevar="rw"
    fi

    cd ..
    return 0    
}

# get an entry from a contrib file, filling the "version" variable
# vwith the version number
#
#   get_contrib_version  contrib_name  file  version
function get_contrib_version(){
    local __resultvar=$3

    # nasty hack: if the name of the "file" is "local_svn", 
    # get the version number  from the local svn checkout of the contribution
    if [[ "$2" == "local_svn" ]]; then
	get_svn_info $1 mode version
	eval $__resultvar="$version"
	return 0
    fi

    # now deal with the version number as if it was an entry in "file" $2
    if [[ -e $2 ]]; then # check if the file actually exists
#      entry=$(grep "^[ \t]*$1[ \t]" $2)  # does not seem to work with tabs
      entry=$(grep "^[[:space:]]*$1[[:space:]]" $2)
      if [ -z "$entry" ]; then
	  eval $__resultvar="[None]"
      else
	  eval $__resultvar="`echo $entry | awk '{print $2}'`"
      fi
    else # file does not exist
      eval $__resultvar="[None]"
    fi    
}

# get a yes/no answer
# returns 0 for n/N/no
#         1 for y/Y/yes
function get_yesno_answer(){
    while true; do
	echo -ne "$1 [y/n] "
	if [[ -z "$2" ]]; then
	    read answer
	else
	    answer="$2"
	    # TODO: add a test that the answer is a valid one
	    echo "$2"
	fi
	case $answer in
	    y|Y|yes) return 1; break ;;
	    n|N|no)  return 0; break ;; 
	esac
    done
}

# check if the local svn has pending modifications
# check_pending_modifications contrib
function check_pending_modifications(){
    cd $1
    result=$(svn status | grep -v "^?")
    if [[ ! -z "$result" ]]; then
	svn status | grep -v "^?"
	cd ..
	return 1
    fi
    cd ..
    return 0
}
