AC_INIT([SISCone], [3.0.5])
AC_CONFIG_SRCDIR([siscone/siscone.cpp])
AM_INIT_AUTOMAKE

dnl (un)comment the following lines if you want to use autoheader
dnl and also (un)comment the call to autoheader in autogen.sh
dnl 
dnl Notes:
dnl  - the first macro needs to place the files in the correct 
dnl    dir in order to have the correct -I options for the remote build
dnl  - see Makefile.am for distcleaning AX_PREFIX_CONFIG_H leftovers
AC_CONFIG_HEADERS(siscone/config_raw.h:config.h.in)
AX_PREFIX_CONFIG_H(siscone/config.h,SISCONE,siscone/config_raw.h)

dnl and also uncomment the call to autoheader in autogen.sh

dnl check autoconf version
AC_PREREQ(2.63)

dnl check basic types
AC_CHECK_TYPE(int)
AC_CHECK_TYPE(long)

dnl set default compilation and link flags
dnl those can be changed at configure time so we don't use AM_CXXFLAGS here
test "x${CXXFLAGS+yes}" = xyes || CXXFLAGS="-O3 -Wall -ffast-math"
dnl CXXFLAGS=" -Wall -O3 -ffast-math "

dnl check useful programs
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

dnl check standard C headers
AC_STDC_HEADERS

dnl set the default destination directory
AC_PREFIX_DEFAULT(/usr/local)

dnl if the debug flag is set, build with -g
dnl default is "yes"
AC_ARG_ENABLE(debug, 
	[  --enable-debug      Turn on debug compiler information],
	[ENABLE_DEBUG_FLAG="$enableval"],
	[ENABLE_DEBUG_FLAG="yes"])
if [[ "x$ENABLE_DEBUG_FLAG" == "xyes" ]] ; then
   CXXFLAGS=${CXXFLAGS}" -g "
fi

dnl uncomment the next line not to build the shared lib by default
dnl AM_DISABLE_SHARED

dnl-----------------
dnl check libraries
dnl ----------------

dnl math lib
AC_CHECK_LIB(m, cos)
dnl already included into LIBS by the previous check
dnl AM_LDFLAGS=" -lm "

dnl if the debug flags are on, check if we can also use 
dnl some profiling tools
dnl COMMENTED: Pass LDFLAGS to configure instead
dnl if [[ "x$ENABLE_DEBUG_FLAG" == "xyes" ]] ; then
dnl    AC_CHECK_LIB(profiler, google_initializer_module_profiler)
dnl    AC_CHECK_LIB(pthread, pthread_create)
dnl    AC_CHECK_LIB(tcmalloc, malloc)
dnl fi

dnl enable libtool
AC_PROG_LIBTOOL

dnl chceck if auto_ptr is deprecated
dnl if yes, see if unique_ptr is supported
AC_MSG_CHECKING([[if $CXX $CXXFLAGS considers atd::auto_ptr as deprecated]])
AC_LANG_PUSH(C++)
save_CXXFLAGS="$CXXFLAGS"
CXXFLAGS="$CXXFLAGS -Werror=deprecated-declarations"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <memory>]],[[int *a= new int(1); std::auto_ptr<int> b; b.reset(a);]])],
    [ac_compilation_deprecated="no"],[ac_compilation_deprecated="yes"])
AC_MSG_RESULT([$ac_compilation_deprecated])
if [[ "$ac_compilation_deprecated" == "yes" ]] ; then
   AC_MSG_CHECKING([[if $CXX $CXXFLAGS supports atd::unique_ptr]])
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <memory>]],[[int *a= new int(1); std::unique_ptr<int> b; b.reset(a);]])],
    [ac_supports_unique_ptr="yes"],[as_supports_unique_ptr="no"])
    if [[ "$ac_supports_unique_ptr" == "yes" ]] ; then
       AC_DEFINE(USES_UNIQUE_PTR_AS_AUTO_PTR, [], [use unique_ptr instead of auto_ptr])
    fi
    AC_MSG_RESULT([$ac_supports_unique_ptr])
    
fi
CXXFLAGS="$save_CXXFLAGS"
AC_LANG_POP(C++)


AC_OUTPUT( Makefile siscone/Makefile siscone/spherical/Makefile examples/Makefile examples/events/Makefile )
