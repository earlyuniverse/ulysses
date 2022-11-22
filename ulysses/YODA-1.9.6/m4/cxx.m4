## Produce more comprehensible error messages if possible. 
## Must be called BEFORE AC_LANG(C++) or it will have no effect.
## See: http://www.bdsoft.com/tools/stlfilt.html
#AC_CEDAR_CXXFILTER()
#----------------------------------------
AC_DEFUN([AC_CEDAR_CXXFILTER], [
  if test "x$CXX" = "x"; then 
    AC_PATH_PROG(GFILT, gfilt, $CXX, $PATH:$HOME/bin:$HOME/local/bin)
    if test "x$GFILT" != "x"; then CXX="$GFILT -banner:N"; fi
  fi
])


## Try to find an "output colorizing" variant on the compiler.
## Haven't decided yet how this should interact with the normal
## AC_PROG_CXX etc.
#AC_CEDAR_CXXCOLOR()
#----------------------------------------
AC_DEFUN([AC_CEDAR_CXXCOLOR], [
  AC_PATH_PROG(COLORCXX, $CXX-color, $CXX, $PATH:$HOME/bin:$HOME/local/bin)
  CXX=$COLORCXX
])


## Determine whether a compiler flag is accepted
#AC_CEDAR_CHECKCXXFLAG(flag, action-if-true, action-if-false)
AC_DEFUN([AC_CEDAR_CHECKCXXFLAG], [
  AC_LANG_PUSH(C++)
  AC_MSG_CHECKING([if the $CXX compiler accepts the $1 flag])
  AC_LANG_CONFTEST([AC_LANG_PROGRAM([],[return 0;])])
  flag_ok=no
  #$CXX $1 conftest.cpp >&5 2>/dev/null && flag_ok=yes
  stat_string=`$CXX $1 conftest.cpp -o yoda-configure.tmp 2>&1 1>&5` ; test -z "$stat_string" && flag_ok=yes
  rm -f yoda-configure.tmp
  AC_MSG_RESULT([$flag_ok])
  if test x$flag_ok = xyes; then 
    true
    $2
  else
    true
    $3
  fi
  AC_LANG_POP(C++)
])
