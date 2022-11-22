# Convert argument to an absolute path
AC_DEFUN([AC_CEDAR_ABSPATH], [
  inpath=${$1}
  abspath=""
  if test -d "$inpath"; then
    abspath=`cd $inpath && pwd`
  elif test -f "$inpath"; then
    dirpart=`dirname $inpath`
    localpart=`basename $inpath`
    abspath="`cd $dirpart && pwd`/$localpart"
  fi
  if test -z "$abspath"; then
    AC_MSG_ERROR(Could not make absolute path $1 from ${$1})
  fi
  $1=$abspath
])


# A version of AC_SUBST which ensures that paths are absolute
AC_DEFUN([AC_CEDAR_PATH_SUBST], [
  AC_CEDAR_ABSPATH($1)
  AC_SUBST($1)
])


# Declares an external package to be used. This macro will declare appropriate 
# command line switches for ./configure and corresponding environment variables
# to define the package's prefix path, library path and include path. The
# configure command line option switches override the corresponding environment
# variables. Various permutations of capitalisation, version number appendage
# and include/lib path will be attempted. Omit the version number if it isn't
# relevant. The third and fourth options will be executed respectively if the
# package is or isn't found: this is an ideal place to place ERROR/WARNING
# messages as appropriate for mandatory/optional packages and to run tests of
# the library if it appears to be there.
#
#AC_CEDAR_LIBRARYANDHEADERS(PrettyName, ReleaseNumber, action-if-true, action-if-false)
#----------------------------------------
AC_DEFUN([AC_CEDAR_LIBRARYANDHEADERS], [
  ## Define a bunch of case permutations
  m4_define([cedar_PkgName], [$1])dnl
  m4_define([cedar_PKGNAME], [translit([translit([$1], [a-z], [A-Z])], [.])])dnl
  m4_define([cedar_pkgname], [translit([translit([$1], [A-Z], [a-z])], [.])])dnl
  m4_define([cedar_SAFEPKGNAME], [translit(cedar_PKGNAME, [-], [_])])dnl
  m4_define([cedar_safepkgname], [translit(cedar_pkgname, [-], [_])])dnl
  m4_define([cedar_pkgversion], [$2])dnl

  pkgpath=""
  pkglib=no; pkginc=no; pkggood=no

  ## Don't know why this isn't working by default:
  test x${prefix} = xNONE && prefix=${ac_default_prefix}

  ## Environment variables for specifying paths
  AC_ARG_VAR(@&t@cedar_SAFEPKGNAME@&t@PATH,
    path to cedar_PkgName @<:@$prefix and various standard locations@:>@)
  pkgpath=${cedar_SAFEPKGNAME@&t@PATH}

  ## "configure" option switches for specifying paths
  AC_ARG_WITH([cedar_safepkgname],
              [AC_HELP_STRING(--with-@&t@cedar_safepkgname, 
                path to cedar_PkgName @<:@$prefix and various standard locations@:>@)],
              [pkgpath=$with_@&t@cedar_safepkgname], [])
  dnl echo "DEBUG: withval=$withval, with_@&t@cedar_safepkgname=$with_@&t@cedar_safepkgname -> pkgpath=$pkgpath"
  if test "$pkgpath"; then cedar_SAFEPKGNAME@&t@PATH="$pkgpath"; fi

  ## Has this lib been disabled?
  #echo "DEBUG: pkgpath = $pkgpath"
  if test x$pkgpath = xno; then 
    AC_MSG_NOTICE(Not building against cedar_PkgName)
    AM_CONDITIONAL(WITH_@&t@cedar_SAFEPKGNAME@&t@INC, false)
    AM_CONDITIONAL(WITH_@&t@cedar_SAFEPKGNAME@&t@LIB, false)
    AM_CONDITIONAL(WITHOUT_@&t@cedar_SAFEPKGNAME@&t@INC, true)
    AM_CONDITIONAL(WITHOUT_@&t@cedar_SAFEPKGNAME@&t@LIB, true)
    $4
  else
    ## Check library and header
    AC_CEDAR_LIBRARY(cedar_PkgName, cedar_pkgversion, , pkglib=no)
    AC_CEDAR_HEADERS(cedar_PkgName, cedar_pkgversion, , pkginc=no)

    ## Execute pass/fail shell code
    if test "x$pkglib" = "xyes" && test "x$pkginc" = "xyes"; then
      #AC_MSG_NOTICE([cedar_PkgName paths verified])
      pkggood="yes"
      $3
    else 
      pkggood="no"
      $4
    fi
  fi

  ## Export variables to automake
  have_@&t@cedar_safepkgname=no
  test x$pkggood != xno && have_@&t@cedar_safepkgname=yes
  AM_CONDITIONAL(WITH_@&t@cedar_SAFEPKGNAME, [test x$pkggood != xno])
  AM_CONDITIONAL(WITHOUT_@&t@cedar_SAFEPKGNAME, [test x$pkggood = xno])
])


#AC_CEDAR_HEADERS(PrettyName, ReleaseNumber, action-if-true, action-if-false)
AC_DEFUN([AC_CEDAR_HEADERS], [
  ## Define a bunch of case permutations
  m4_define([cedar_PkgName], [$1])dnl
  m4_define([cedar_PKGNAME], [translit([translit([$1], [a-z], [A-Z])], [.])])dnl
  m4_define([cedar_pkgname], [translit([translit([$1], [A-Z], [a-z])], [.])])dnl
  m4_define([cedar_SAFEPKGNAME], [translit(cedar_PKGNAME, [-], [_])])dnl
  m4_define([cedar_safepkgname], [translit(cedar_pkgname, [-], [_])])dnl
  m4_define([cedar_IncName], [cedar_PkgName])dnl
  m4_define([cedar_INCNAME], [cedar_PKGNAME])dnl
  m4_define([cedar_incname], [cedar_pkgname])dnl
  m4_define([cedar_IncName1], [translit(cedar_PkgName, [-], [_])])dnl
  m4_define([cedar_INCNAME1], [translit(cedar_PKGNAME, [-], [_])])dnl
  m4_define([cedar_incname1], [translit(cedar_pkgname, [-], [_])])dnl
  m4_define([cedar_IncName2], [translit(cedar_PkgName, [-], [])])dnl
  m4_define([cedar_INCNAME2], [translit(cedar_PKGNAME, [-], [])])dnl
  m4_define([cedar_incname2], [translit(cedar_pkgname, [-], [])])dnl
  m4_define([cedar_pkgversion], [$2])dnl

  ## We have a set of user-set variables:
  pkgpath=""; pkgincpath=""
  ## Also need some status variables:
  pkginc=no

  ## Don't know why this isn't working by default:
  test x${prefix} = xNONE && prefix=${ac_default_prefix}

  ## Environment variables for specifying paths
  AC_ARG_VAR(@&t@cedar_SAFEPKGNAME@&t@PATH,
    path to cedar_PkgName @<:@$prefix and various standard locations@:>@)
  AC_ARG_VAR(@&t@cedar_SAFEPKGNAME@&t@INCPATH,
    path to the directory containing the cedar_PkgName header files @<:@cedar_SAFEPKGNAME@&t@PATH/include@:>@)

  ## "configure" option switches for specifying paths
  pkgpath=${cedar_SAFEPKGNAME@&t@PATH}
  AC_ARG_WITH(cedar_safepkgname,
              AC_HELP_STRING(--with-@&t@cedar_safepkgname@&t@, 
                path to cedar_PkgName @<:@$prefix and various standard locations@:>@),
              [pkgpath=$with_@&t@cedar_safepkgname], [])
  dnl echo "DEBUG: withval=$withval, with_@&t@cedar_safepkgname=$with_@&t@cedar_safepkgname -> pkgpath=$pkgpath"
  if test "$pkgpath"; then cedar_SAFEPKGNAME@&t@PATH="$pkgpath"; fi

  pkgincpath=${cedar_SAFEPKGNAME@&t@INCPATH}
  AC_ARG_WITH(cedar_safepkgname@&t@-incpath,
              AC_HELP_STRING(--with-@&t@cedar_safepkgname@&t@-incpath, 
                path to directory containing cedar_PkgName headers @<:@cedar_SAFEPKGNAME@&t@PATH/include@:>@),
              [pkgincpath=$with_@&t@cedar_safepkgname@&t@_incpath], [])

  ## Has this header been disabled?
  if test x$pkgpath = xno; then 
    AC_MSG_NOTICE(Not building against cedar_PkgName)
    $4
  else
    ## Base paths
    pkgbases="$prefix $ac_default_prefix /usr /"
    if test "$pkgpath"; then pkgbases="$pkgpath"; fi

    ## Look for include files: first build the search list...
    incpaths=""
    if test "$pkgincpath"; then 
      incpath=`echo $pkgincpath | sed -e 's://*:/:g' -e 's:/$::'`
      incpaths="$incpath"
    else
      for base in $pkgbases; do
        incpath=`echo "$base/include" | sed -e 's://*:/:g' -e 's:/$::'`
        incpaths="$incpaths $incpath"
      done
    fi
    
    if test "x$CEDAR_M4_DEBUG" != "x"; then
      echo "DEBUG: inc paths = $incpaths"
    fi


    ## Build package names
    incnames="cedar_IncName cedar_INCNAME cedar_incname"
    if test "cedar_IncName" != "cedar_IncName1"; then
      incnames="$incnames cedar_IncName1 cedar_INCNAME1 cedar_incname1"
    fi
    if test "cedar_IncName" != "cedar_IncName2"; then
      incnames="$incnames cedar_IncName2 cedar_INCNAME2 cedar_incname2"
    fi

    ## .. and then do the search:
    for incpath in $incpaths; do
      for incname in $incnames; do
        fullincpath="$incpath/$incname"
        if test "x$CEDAR_M4_DEBUG" != "x"; then
          echo "Testing cedar_PkgName inc path: $fullincpath"
        fi
        if test -d $fullincpath; then
          pkginc=yes
          break
        else 
          pkginc=no;
        fi
        if test x$pkginc != xno; then break; fi
      done
      if test x$pkginc != xno; then break; fi
    done

    if test x$pkginc != xno; then
      cedar_SAFEPKGNAME@&t@INCPATH="$incpath"
      AC_CEDAR_ABSPATH(cedar_SAFEPKGNAME@&t@INCPATH)
      cedar_SAFEPKGNAME@&t@INCNAME="$incname"
      cedar_SAFEPKGNAME@&t@CPPFLAGS="-I$cedar_SAFEPKGNAME@&t@INCPATH"
      #echo cedar_SAFEPKGNAME@&t@INCPATH : $cedar_SAFEPKGNAME@&t@INCPATH
      AC_MSG_NOTICE([Found cedar_PkgName header directory at $incpath])
    else
      ## Last resort --- only tried if $pkgpath was specified
      if test x$pkgpath != x; then
        incpath="$pkgpath/include"
        if test -d "$incpath"; then
          cedar_SAFEPKGNAME@&t@INCPATH=`echo $incpath | sed -e s:'/$':'':`
          AC_CEDAR_ABSPATH(cedar_SAFEPKGNAME@&t@INCPATH)
          cedar_SAFEPKGNAME@&t@INCNAME=""
          cedar_SAFEPKGNAME@&t@CPPFLAGS="-I$cedar_SAFEPKGNAME@&t@INCPATH"
          pkginc=yes
          #AC_MSG_NOTICE([Found cedar_PkgName header directory at $incpath])
        fi
      else
        AC_MSG_WARN(cedar_PkgName header directory was not found)
      fi
    fi

    ## Execute pass/fail shell code
    if test "x$pkginc" = "xyes"; then
      true
      $3
    else 
      true
      $4
    fi
  fi

  ## Export variables to automake
  AC_SUBST(cedar_SAFEPKGNAME@&t@INCPATH)
  AC_SUBST(cedar_SAFEPKGNAME@&t@INCNAME)
  AC_SUBST(cedar_SAFEPKGNAME@&t@CPPFLAGS)
  AM_CONDITIONAL(WITH_@&t@cedar_SAFEPKGNAME@&t@INC, [test x$pkgginc != xno])
  AM_CONDITIONAL(WITHOUT_@&t@cedar_SAFEPKGNAME@&t@INC, [test x$pkgginc = xno])
])


#AC_CEDAR_LIBRARY(PrettyName, ReleaseNumber, action-if-true, action-if-false)
AC_DEFUN([AC_CEDAR_LIBRARY], [
  ## Define a bunch of case permutations
  m4_define([cedar_PkgName], [$1])dnl
  m4_define([cedar_PKGNAME], [translit([translit([$1], [a-z], [A-Z])], [.])])dnl
  m4_define([cedar_pkgname], [translit([translit([$1], [A-Z], [a-z])], [.])])dnl
  m4_define([cedar_SAFEPKGNAME], [translit(cedar_PKGNAME, [-], [_])])dnl
  m4_define([cedar_safepkgname], [translit(cedar_pkgname, [-], [_])])dnl
  m4_define([cedar_LibName], [cedar_PkgName])dnl
  m4_define([cedar_LIBNAME], [cedar_PKGNAME])dnl
  m4_define([cedar_libname], [cedar_pkgname])dnl
  m4_define([cedar_LibName1], [translit(cedar_PkgName, [-], [_])])dnl
  m4_define([cedar_LIBNAME1], [translit(cedar_PKGNAME, [-], [_])])dnl
  m4_define([cedar_libname1], [translit(cedar_pkgname, [-], [_])])dnl
  m4_define([cedar_LibName2], [translit(cedar_PkgName, [-], [])])dnl
  m4_define([cedar_LIBNAME2], [translit(cedar_PKGNAME, [-], [])])dnl
  m4_define([cedar_libname2], [translit(cedar_pkgname, [-], [])])dnl
  m4_define([cedar_libversion], [$2])dnl

  ## We have a set of user-set variables:
  pkgpath=""; pkglibpath=""; pkglibname=""
  ## Also need some status variables:
  pkglib=no

  ## Don't know why this isn't working by default:
  test x${prefix} = xNONE && prefix=${ac_default_prefix}

  ## Environment variables for specifying paths
  AC_ARG_VAR(@&t@cedar_SAFEPKGNAME@&t@PATH,
    path to cedar_PkgName @<:@$prefix and various standard locations@:>@)
  AC_ARG_VAR(@&t@cedar_SAFEPKGNAME@&t@LIBPATH,
    path to the directory containing the cedar_PkgName library @<:@cedar_SAFEPKGNAME@&t@PATH/lib or cedar_SAFEPKGNAME@&t@PATH/lib/cedar_PkgName@:>@)
  AC_ARG_VAR(@&t@cedar_SAFEPKGNAME@&t@LIBNAME,
    name to be used when linking the cedar_PkgName library @<:@cedar_PkgName@:>@)
  pkgpath=${cedar_SAFEPKGNAME@&t@PATH}

  ## "configure" option switches for specifying paths
  AC_ARG_WITH([cedar_safepkgname],
              [AC_HELP_STRING(--with-@&t@cedar_safepkgname, 
                path to cedar_PkgName @<:@$prefix and various standard locations@:>@)],
              [pkgpath=$with_@&t@cedar_safepkgname], [])

  dnl echo "DEBUG: withval=$withval, with_@&t@cedar_safepkgname=$with_@&t@cedar_safepkgname -> pkgpath=$pkgpath"
  if test "$pkgpath"; then cedar_SAFEPKGNAME@&t@PATH="$pkgpath"; fi
  pkglibpath=${cedar_SAFEPKGNAME@&t@LIBPATH}
  pkglibname=${cedar_SAFEPKGNAME@&t@LIBFLAG}

  AC_ARG_WITH(cedar_safepkgname@&t@-libpath,
              AC_HELP_STRING(--with-@&t@cedar_safepkgname@&t@-libpath, 
                path to directory containing cedar_PkgName library @<:@cedar_SAFEPKGNAME@&t@PATH/lib or cedar_SAFEPKGNAME@&t@PATH/lib/cedar_PkgName@:>@),
              [pkglibpath=$with_@&t@cedar_safepkgname@&t@_libpath], [])
  AC_ARG_WITH(cedar_safepkgname@&t@-libname,
              AC_HELP_STRING(--with-@&t@cedar_safepkgname@&t@-libname,
                name to be used when linking the cedar_PkgName library @<:@cedar_PkgName@:>@),
              [pkglibname=$with_@&t@cedar_safepkgname@&t@_libname], [])

  ## Has this lib been disabled?
  if test x$pkgpath = xno; then 
    AC_MSG_NOTICE(Not building against cedar_PkgName)
    $4
  else
    ## Base paths
    pkgbases="$prefix $ac_default_prefix /usr /"
    if test "$pkgpath"; then pkgbases="$pkgpath"; fi
    if test "x$CEDAR_M4_DEBUG" != "x"; then
      echo "DEBUG: $pkgpath -> $pkgbases"
    fi

    ## Build a list of library search locations, unless specified
    libdirnames="lib"

    ## Test for 64-bit mode and add lib64 as first choice 
    ## library dir name if appropriate
    if test -z "$UNAME"; then 
      AC_PATH_PROG(UNAME, [uname], [no])
    fi
    if test x$UNAME != xno; then
      if test -n `$UNAME -m | grep 64`; then
        libdirnames="lib64 $libdirnames"
      fi
    fi

    libpaths=""
    if test "$pkglibpath"; then 
      libpath=`echo $pkglibpath | sed -e 's://*:/:g' -e 's:/$::'`
      libpaths="$libpath"
    else
      ## Outer loop over lib / lib64 part
      for libdirname in $libdirnames; do 
        ## Inner loop over base path
        for base in $pkgbases; do
          libpath=`echo "$base/$libdirname" | sed -e 's://*:/:g' -e 's:/$::'`
          libpaths="$libpaths $libpath"
        done
      done
      libpaths="$libpaths ./src"
    fi

    ## Use case permuatations on the package name (mixed, all-upper and 
    ## all-lower), plus punctuation replacements.
    libnames="cedar_LibName cedar_LIBNAME cedar_libname"
    if test "cedar_LibName" != "cedar_LibName1"; then
      libnames="$libnames cedar_LibName1 cedar_LIBNAME1 cedar_libname1"
    fi
    if test "cedar_LibName" != "cedar_LibName2"; then
      libnames="$libnames cedar_LibName2 cedar_LIBNAME2 cedar_libname2"
    fi
    if test "$pkglibname"; then libnames=$pkglibname; fi

    ## Define library versions
    libversions=""
    if test x"cedar_libversion" != x; then
      libversions="-cedar_libversion cedar_libversion"
    fi

    ## Look for library with various name permutations  
    for libpath in $libpaths; do
      for libversion in $libversions ""; do
        for libname in $libnames; do
          for libextn in la so dylib dll a; do
            testpath="${libpath}/lib${libname}${libversion}.${libextn}"
            testpath=`echo $testpath | sed -e 's://*:/:g' -e 's:/$::'`
            if test "x$CEDAR_M4_DEBUG" != "x"; then
              echo "DEBUG: Testing $testpath"
            fi
            if test -e $testpath; then 
              pkglib=yes
              break
            else
              pkglib=no
            fi
          done
          if test x$pkglib != xno; then break; fi
        done
        if test x$pkglib != xno; then break; fi
      done
      if test x$pkglib != xno; then break; fi
    done

    ## Announce success/failure and set variables
    if test x$pkglib != xno; then
      ## TODO: Can we just use $testpath here?
      libfullpath=`echo ${libpath}/lib${libname}${libversion}.${libextn}`
      libfullpath=`echo $libfullpath | sed -e 's://*:/:g' -e 's:/$::'`
      AC_CEDAR_ABSPATH(libfullpath)
      cedar_SAFEPKGNAME@&t@LIB="$libfullpath"
      cedar_SAFEPKGNAME@&t@LIBPATH=`dirname $libfullpath`
      ## TODO: Can we just use $libname here?
      cedar_SAFEPKGNAME@&t@LIBNAME=`basename $libfullpath | sed -e s/'^lib'// -e s/'\.@<:@a-zA-Z@:>@*$'//`
      cedar_SAFEPKGNAME@&t@LDFLAGS="-L$cedar_SAFEPKGNAME@&t@LIBPATH"
      cedar_SAFEPKGNAME@&t@LDLIBS="-l$cedar_SAFEPKGNAME@&t@LIBNAME"
      AC_MSG_NOTICE([Found cedar_PkgName library at ${libfullpath}])
    else
      AC_MSG_WARN(cedar_PkgName library was not found)
    fi

    ## Execute pass/fail shell code
    if test "x$pkglib" = "xyes"; then
      true
      $3
    else 
      true
      $4
    fi
  fi

  ## Export variables to automake
  AC_SUBST(cedar_SAFEPKGNAME@&t@LIB)
  AC_SUBST(cedar_SAFEPKGNAME@&t@LIBPATH)
  AC_SUBST(cedar_SAFEPKGNAME@&t@LIBNAME)
  AC_SUBST(cedar_SAFEPKGNAME@&t@LDFLAGS)
  AC_SUBST(cedar_SAFEPKGNAME@&t@LDLIBS)
  AM_CONDITIONAL(WITH_@&t@cedar_SAFEPKGNAME@&t@LIB, [test x$pkglib != xno])
  AM_CONDITIONAL(WITHOUT_@&t@cedar_SAFEPKGNAME@&t@LIB, [test x$pkglib = xno])
])


# Create a variable called EMPTY with appropriate content
#AC_EMPTY_SUBST
AC_DEFUN([AC_EMPTY_SUBST],
[EMPTY=""
AC_SUBST(EMPTY)
])
