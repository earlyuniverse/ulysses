#AC_CEDAR_OSX()
#----------------------------------------
AC_DEFUN([AC_CEDAR_OSX], [
  LIBPATHVARNAME="LD_LIBRARY_PATH"
  AC_CHECK_TOOL(SWVERS, sw_vers)
  if test x$SWVERS != x; then
    PROD_NAME=$($SWVERS -productName | cut -f 2 -d:)
  fi
  AM_CONDITIONAL(WITH_OSX, [test "$PROD_NAME" = "Mac OS X" -o "$PROD_NAME" = "macOS"])
  if test "$PROD_NAME" = "Mac OS X" -o "$PROD_NAME" = "macOS"; then
    #MACOSX_DEPLOYMENT_TARGET=$($SWVERS -productVersion | cut -f 1,2 -d.)
    #AC_MSG_NOTICE([MACOSX_DEPLOYMENT_TARGET = $MACOSX_DEPLOYMENT_TARGET])
    AM_CXXFLAGS="$AM_CXXFLAGS -Dunix"
    LIBPATHVARNAME="DYLD_LIBRARY_PATH"
  fi
  AC_SUBST(LIBPATHVARNAME)
])
