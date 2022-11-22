AC_INIT([FastJet],[3.4.0])
AC_CONFIG_SRCDIR([src/JetDefinition.cc])
AM_INIT_AUTOMAKE(tar-ustar)

dnl Make sure users can't install in build directory
if test "$prefix" = "$PWD"; then
  AC_MSG_ERROR([Installation into the build directory is not supported: use a different --prefix argument])
fi

dnl (un_comment the following lines if you want to use autoheader
dnl and also (un_comment the call to autoheader in autogen.sh
dnl 
dnl Notes:
dnl  - the first macro needs to place the files in the correct 
dnl    dir in order to have the correct -I options for the remote build
dnl  - see Makefile.am for distcleaning AX_PREFIX_CONFIG_H leftovers
AC_CONFIG_HEADERS(include/fastjet/config_raw.h:config.h.in)
AX_PREFIX_CONFIG_H(include/fastjet/config_auto.h,FASTJET,include/fastjet/config_raw.h)

dnl check autoconf version
dnl
dnl 2011-07-07 tested explicitly to be working with
dnl   autoconf-2.59 (lxplus)
dnl   autoconf-2.61 (karnak)
dnl   autoconf-2.65 (ipht)
dnl   autoconf-2.68 (virtual&lpthe)
dnl
dnl 2013-03-30, updated and tested it works with
dnl   autoconf-2.63 + automake 1.11 (tycho SLC6)
dnl   autoconf-2.69 + automake 1.13 (gavin's laptop with macports)
dnl
dnl   NB: it is broken on lxplus (SLC5, autoconf 2.59 + automake 1.9), 
dnl       but lxplus will soon migrate to SLC6. For now 
dnl       autogen can be run on lxplus6, links to automake-1.11 
dnl       files replaced with copies and then used on lxplus5
AC_PREREQ(2.63)

dnl load our M4 macro directory
dnl  - this is recommended by latest versions of libtoolize
dnl  - this is equivalent to adding ACLOCAL_AMFLAGS = -I m4 in Makefile.am
dnl  - the libtoolize manpage says that for the moment it's more 
dnl    portable to keep using ACLOCAL_AMFLAGS.
dnl  - the autoconf manual mentions that since we use 'aclocal -I m4'
dnl    we need ACLOCAL_AMFLAGS anyway
dnl as a consequence, we don't use AC_CONFIG_MACRO_DIR for the moment
dnl Question: is it compatible with older versions of libtool?
dnl   AC_CONFIG_MACRO_DIR([m4])

dnl check basic types
AC_CHECK_TYPE(int)
AC_CHECK_TYPE(long)

dnl-------------------------------------------------------------------
dnl if the user has not set CXXFLAGS, then set our own default
dnl-------------------------------------------------------------------
USE_CUSTOM_FLAGS="yes"
if [[ "x$CXXFLAGS" == "x" ]]; then
   dnl 2014-08-11: it seems that -O3 produces timing results that vary
   dnl quite a lot between different architectures and are sometimes
   dnl worse than -O2. So we'll use -O2 as a new default
   CXXFLAGS="-O2 -Wall"
   USE_CUSTOM_FLAGS="no"
fi

# test "x${CXXFLAGS+yes}" == "xyes" || CXXFLAGS="-O3 -Wall"

dnl-------------------------------------------------------------------
dnl check useful programs
dnl-------------------------------------------------------------------
AC_PROG_CXX

dnl ensure CXX presence
cxx_progname=$(basename $CXX)
cxx_dirname=$(dirname $CXX)
AC_CHECK_PROG([cxx_found],[$cxx_progname],[yes],[no],[$cxx_dirname$PATH_SEPARATOR$PATH])
if test "x$cxx_found" == "xno"; then
   AC_MSG_ERROR([No functional C++ compiler found])
fi

AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_SED
AC_PROG_MAKE_SET

AM_CXXFLAGS=""

dnl check standard C headers
AC_STDC_HEADERS

dnl set the default destination directory
AC_PREFIX_DEFAULT(/usr/local)

dnl-------------------------------------------------------------------
dnl initialize compilation and link flags to set in fastjet-config
dnl-------------------------------------------------------------------
CONFIG_CXXFLAGS=""
CONFIG_RUNPATH_FLAGS=""
CONFIG_LIBS=""
CONFIG_LIBS_PLUGINS=""
CONFIG_LIBS_PLUGINS_STATIC=""

dnl-------------------------------------------------------------------
dnl deal with debug option
dnl-------------------------------------------------------------------
dnl if the debug flag is set, build with -g
dnl default is "yes"
AC_ARG_ENABLE(debug,
	[  --enable-debug          Turn on debug compiler information [default=yes]],
	[ENABLE_DEBUG_FLAG="$enableval"],
	[ENABLE_DEBUG_FLAG="yes"])
if [[ "$ENABLE_DEBUG_FLAG" == "yes" ]] ; then
   AM_CXXFLAGS=${AM_CXXFLAGS}" -g "

   dnl-------------------------------------------------------------------
   dnl check if the compiler supports -Woverloaded-virtual
   dnl-------------------------------------------------------------------
   AC_MSG_CHECKING([if $CXX supports -Woverloaded-virtual flag])
   SAVE_CXXFLAGS="$CXXFLAGS"
   CXXFLAGS="-Woverloaded-virtual"
   AC_LANG_PUSH(C++)
   AC_COMPILE_IFELSE([AC_LANG_SOURCE([int main(void){ return 0;}])],
       [ac_overloaded_virtual_ok=yes],[ac_overloaded_virtual_ok=no])
   AC_LANG_POP(C++)
   CXXFLAGS="$SAVE_CXXFLAGS"
   AC_MSG_RESULT([$ac_overloaded_virtual_ok])
   if [[ "x${ac_overloaded_virtual_ok}" == "xyes" ]]; then
      AM_CXXFLAGS=${AM_CXXFLAGS}" -Woverloaded-virtual "
   fi
fi

dnl check if one wants a high level of warnings
AC_ARG_ENABLE(extra-warnings,
	[  --enable-extra-warnings Turn on additional compiler warnings [default=no]],
	[ENABLE_EXTRA_WARNINGS_FLAG="$enableval"],
	[ENABLE_EXTRA_WARNINGS_FLAG="no"])
if [[ "$ENABLE_EXTRA_WARNINGS_FLAG" == "yes" ]] ; then
   AM_CXXFLAGS=${AM_CXXFLAGS}" -pedantic -Wextra -Wshadow "
fi


dnl-------------------------------------------------------------------
dnl check if a few modern compiler options are availabble

dnl --------- disable deprecated auto_ptr interface ------------------
dnl This requires the --disable-auto-ptr option 
dnl which is currently not the default but will become in a future
dnl release (or might disappear altogether)
AC_ARG_ENABLE(auto-ptr,
	[  --enable-auto-ptr  compile the deprecated parts of the interface using auto-ptr [default=yes]],
	[ENABLE_AUTO_PTR_INTERFACE="$enableval"],
	[ENABLE_AUTO_PTR_INTERFACE="yes"])
if [[ "$ENABLE_AUTO_PTR_INTERFACE" == "yes" ]] ; then
   AC_DEFINE(HAVE_AUTO_PTR_INTERFACE, [], [compile the deprecated parts of the interface using auto-ptr])   
fi

dnl check if one wants to enable c++11 features (required e.g. for thread-safety)
dnl
dnl This defines a FASTJET_HAVE_CXX11_FEATURES in fastjet/config.h
AC_ARG_ENABLE(limited-thread-safety,
	[  --enable-limited-thread-safety   Turn on basic thread-safety features (requires C++11)  [default=no]],
	[ENABLE_LIMITED_THREAD_SAFETY_FLAG="$enableval"],
	[ENABLE_LIMITED_THREAD_SAFETY_FLAG="no"])
if [[ "$ENABLE_LIMITED_THREAD_SAFETY_FLAG" == "yes" ]] ; then
   AC_DEFINE(HAVE_LIMITED_THREAD_SAFETY, [], [defined if limited thread-safety has been enabled])
   AC_MSG_WARN([Thread safety uses C++11 features and may require additional compiler arguments, these have to be specified by the user])
   # one could alternately try to see what compiler flag(s) has to be included
   # although that can be fairly heavy so better be moved in an ac file 
   #
   # e.g. for g++ we could try
   #   AC_MSG_CHECKING([if $CXX supports -std=c++11 flag])
   #   SAVE_CXXFLAGS="$CXXFLAGS"
   #   CXXFLAGS="-std=c++11"
   #   AC_LANG_PUSH(C++)
   #   AC_COMPILE_IFELSE([AC_LANG_SOURCE([int main(void){ return 0;}])],
   #       [ac_cxx_eleven_ok=yes],[ac_cxx_eleven_ok=no])
   #   AC_LANG_POP(C++)
   #   CXXFLAGS="$SAVE_CXXFLAGS"
   #   AC_MSG_RESULT([$ac_cxx_eleven_ok])
   #   if [[ "x${ac_cxx_eleven_ok}" == "xyes" ]]; then
   #      AM_CXXFLAGS=${AM_CXXFLAGS}" -std=c++11 "
   #   fi
fi

AC_ARG_ENABLE(thread-safety,
	[  --enable-thread-safety           Turn on thread-safety features (requires C++11)  [default=no]],
	[ENABLE_THREAD_SAFETY_FLAG="$enableval"],
	[ENABLE_THREAD_SAFETY_FLAG="no"])
if [[ "$ENABLE_THREAD_SAFETY_FLAG" == "yes" ]] ; then
   AC_DEFINE(HAVE_LIMITED_THREAD_SAFETY, [], [defined if limited thread-safety has been enabled])
   AC_DEFINE(HAVE_THREAD_SAFETY,         [], [defined if thread-safety has been enabled])
   AC_MSG_WARN([Thread safety uses C++11 features and may require additional compiler arguments, these have to be specified by the user])
fi

dnl-------------------------------------------------------------------
dnl check if a few modern compiler options are availabble

dnl---------- [[deprecated]] --> HAVE_CXX14_DEPRECATED ---------------
AC_MSG_CHECKING([[if $CXX $CXXFLAGS supports the "[[deprecated]]" keyword]])
AC_LANG_PUSH(C++)
save_CXXFLAGS="$CXXFLAGS"
CXXFLAGS="$CXXFLAGS -Werror=attributes"
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[class A{public:  [[deprecated]] void test() {}}; int main(void){ A a; a.test(); return 0;}]])],
    [ac_compiler_has_requested_support="yes"],[ac_compiler_has_requested_support="no"])
CXXFLAGS="$save_CXXFLAGS"
AC_LANG_POP(C++)
if [[ "$ac_compiler_has_requested_support" == "yes" ]] ; then
   AC_DEFINE(HAVE_CXX14_DEPRECATED, [], [compiler supports c++14 deprecated keyword])
fi
AC_MSG_RESULT([$ac_compiler_has_requested_support])

dnl------ _attribute__((deprecated)) --> HAVE_GNUCXX_DEPRECATED ------
AC_MSG_CHECKING([[if $CXX $CXXFLAGS supports the "__attribute__((deprecated))" keyword]])
AC_LANG_PUSH(C++)
save_CXXFLAGS="$CXXFLAGS"
CXXFLAGS="$CXXFLAGS -Werror=attributes"
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[class A{public: void test() __attribute__((deprecated)) {}};int main(void){ A a; a.test(); return 0;}]])],
    [ac_compiler_has_requested_support="yes"],[ac_compiler_has_requested_support="no"])
CXXFLAGS="$save_CXXFLAGS"
AC_LANG_POP(C++)
if [[ "$ac_compiler_has_requested_support" == "yes" ]] ; then
   AC_DEFINE(HAVE_GNUCXX_DEPRECATED, [], [compiler supports GNU c++ deprecated attribute])
fi
AC_MSG_RESULT([$ac_compiler_has_requested_support])

dnl------------------- override --> HAVE_OVERRIDE ---------------------
AC_MSG_CHECKING([[if $CXX $CXXFLAGS supports the "override" keyword]])
AC_LANG_PUSH(C++)
save_CXXFLAGS="$CXXFLAGS"
CXXFLAGS="$CXXFLAGS -Werror"
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[class A{public: virtual void test(){}}; class B : public A{public: virtual void test() override{}}; int main(void){ B b; b.test(); return 0;}]])],
    [ac_compiler_has_requested_support="yes"],[ac_compiler_has_requested_support="no"])
CXXFLAGS="$save_CXXFLAGS"
AC_LANG_POP(C++)
if [[ "$ac_compiler_has_requested_support" == "yes" ]] ; then
   AC_DEFINE(HAVE_OVERRIDE, [], [compiler supports the "override" keyword])
fi
AC_MSG_RESULT([$ac_compiler_has_requested_support])

dnl---------- [[deprecated]] --> HAVE_CXX14_DEPRECATED ---------------
AC_MSG_CHECKING([[if $CXX ${AM_CXXFLAGS} $CXXFLAGS supports the "[[deprecated]]" keyword]])
AC_LANG_PUSH(C++)
save_CXXFLAGS="$CXXFLAGS"
CXXFLAGS="${AM_CXXFLAGS} $CXXFLAGS -Werror -Wno-error=deprecated-declarations"
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[class A{public:  [[deprecated]] void test() {}}; int main(void){ A a; a.test(); return 0;}]])],
    [ac_compiler_has_requested_support="yes"],[ac_compiler_has_requested_support="no"])
CXXFLAGS="$save_CXXFLAGS"
AC_LANG_POP(C++)
if [[ "$ac_compiler_has_requested_support" == "yes" ]] ; then
   AC_DEFINE(HAVE_CXX14_DEPRECATED, [], [compiler supports c++14 deprecated keyword])
fi
AC_MSG_RESULT([$ac_compiler_has_requested_support])

dnl------ _attribute__((deprecated)) --> HAVE_GNUCXX_DEPRECATED ------
AC_MSG_CHECKING([[if $CXX ${AM_CXXFLAGS} $CXXFLAGS supports the "__attribute__((deprecated))" keyword]])
AC_LANG_PUSH(C++)
save_CXXFLAGS="$CXXFLAGS"
CXXFLAGS="${AM_CXXFLAGS} $CXXFLAGS -Werror -Wno-error=deprecated-declarations"
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[class A{public: __attribute__((__deprecated__)) void test() {}};int main(void){ A a; a.test(); return 0;}]])],
    [ac_compiler_has_requested_support="yes"],[ac_compiler_has_requested_support="no"])
CXXFLAGS="$save_CXXFLAGS"
AC_LANG_POP(C++)
if [[ "$ac_compiler_has_requested_support" == "yes" ]] ; then
   AC_DEFINE(HAVE_GNUCXX_DEPRECATED, [], [compiler supports GNU c++ deprecated attribute])
fi
AC_MSG_RESULT([$ac_compiler_has_requested_support])

dnl------------------- override --> HAVE_OVERRIDE ---------------------
AC_MSG_CHECKING([[if $CXX ${AM_CXXFLAGS} $CXXFLAGS supports the "override" keyword]])
AC_LANG_PUSH(C++)
save_CXXFLAGS="$CXXFLAGS"
CXXFLAGS="${AM_CXXFLAGS} $CXXFLAGS -Werror"
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[class A{public: virtual void test(){}}; class B : public A{public: virtual void test() override{}}; int main(void){ B b; b.test(); return 0;}]])],
    [ac_compiler_has_requested_support="yes"],[ac_compiler_has_requested_support="no"])
CXXFLAGS="$save_CXXFLAGS"
AC_LANG_POP(C++)
if [[ "$ac_compiler_has_requested_support" == "yes" ]] ; then
   AC_DEFINE(HAVE_OVERRIDE, [], [compiler supports the "override" keyword])
fi
AC_MSG_RESULT([$ac_compiler_has_requested_support])

dnl--- explicit (for operators in class) --> HAVE_EXPLICIT_FOR_OPERATORAS ---
AC_MSG_CHECKING([[if $CXX ${AM_CXXFLAGS} $CXXFLAGS supports the "explicit" keyword for operators]])
AC_LANG_PUSH(C++)
save_CXXFLAGS="$CXXFLAGS"
CXXFLAGS="${AM_CXXFLAGS} $CXXFLAGS -Werror"
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[class A{public: explicit operator bool(){return true;}}; int main(void){ A a; bool b=bool(a); if (b) b=false; return 0;}]])],
    [ac_compiler_has_requested_support="yes"],[ac_compiler_has_requested_support="no"])
CXXFLAGS="$save_CXXFLAGS"
AC_LANG_POP(C++)
if [[ "$ac_compiler_has_requested_support" == "yes" ]] ; then
   AC_DEFINE(HAVE_EXPLICIT_FOR_OPERATORS, [], [compiler supports the "explicit" keyword for operators])
fi
AC_MSG_RESULT([$ac_compiler_has_requested_support])



dnl-------------------------------------------------------------------
dnl deal with the choice of shared/static libs
dnl-------------------------------------------------------------------
dnl build the shared lib by default
dnl even though they're enabled by default, AM_ENABLE_SHARED has to
dnl be called so that ${enable_shared} is set to 'yes' instead of ''
AM_ENABLE_SHARED
dnl
dnl comment the previous line and uncomment this one if you want shared 
dnl libs to be disabled by default
dnl AM_DISABLE_SHARED
HAS_SHARED=${enable_shared}

dnl-------------------------------------------------------------------
dnl check headers
dnl-------------------------------------------------------------------
dnl execinfo.h (use in Error.cc to have the backtrace)
AC_CHECK_HEADERS(execinfo.h)

dnl-------------------------------------------------------------------
dnl check libraries
dnl-------------------------------------------------------------------
dnl math lib
AC_CHECK_LIB(m, cos)
CONFIG_CXXFLAGS=""
CONFIG_LIBS=${CONFIG_LIBS}" -lm "

dnl-------------------------------------------------------------------
dnl deal with CGAL
dnl-------------------------------------------------------------------
CGAL_CPPFLAGS=""
CGAL_LIBS=""

dnl first check if the user wants to use CGAL
dnl default value is "no"
dnl Note that CGAL compilation flags are considered as required
dnl so we do store them in "AM_CXXFLAGS"
AC_ARG_ENABLE(cgal,
       [  --enable-cgal           enables link with the CGAL library [default=no]],
       [ENABLE_CGAL_FLAG="$enableval"],
       [ENABLE_CGAL_FLAG="no"])
dnl if asked for, check if the CGAL library can be used
if [[ "x${ENABLE_CGAL_FLAG}" == "xyes" ]] ; then
   dnl if the search for CGAL has been explicitly asked and not found, return an error
   dnl ACX_CHECK_CGAL([],[ENABLE_CGAL_FLAG="no"])
   ACX_CHECK_CGAL([],AC_MSG_ERROR([CGAL explicitly requested and not found or not functional]))
fi

dnl update compilation/link flags for the CGAL lib
if [[ "x$ENABLE_CGAL_FLAG" == "xyes" ]] ; then
   AM_CXXFLAGS=${AM_CXXFLAGS}" $CGAL_CPPFLAGS "
   CONFIG_LIBS=" $CGAL_LIBS "${CONFIG_LIBS}
else
   AM_CXXFLAGS=${AM_CXXFLAGS}" -DDROP_CGAL "
dnl   CONFIG_CXXFLAGS=${CONFIG_CXXFLAGS}" -DDROP_CGAL "
fi
AC_SUBST(CGAL_CPPFLAGS)
AC_SUBST(CGAL_LIBS)

dnl-------------------------------------------------------------------
dnl check if the linker supports --enable-new-dtags (rpath-related)
dnl-------------------------------------------------------------------
AC_MSG_CHECKING([if $CXX supports -Wl,--enable-new-dtags flag])
SAVE_CXXFLAGS="$CXXFLAGS"
CXXFLAGS="-Werror -Wl,--enable-new-dtags"
AC_LANG_PUSH(C++)
AC_LINK_IFELSE([AC_LANG_SOURCE([int main(void){ return 0;}])],
   [ac_enable_new_dtags_ok=yes],[ac_enable_new_dtags_ok=no])
AC_LANG_POP(C++)
CXXFLAGS="$SAVE_CXXFLAGS"
AC_MSG_RESULT([$ac_enable_new_dtags_ok])
HAS_RUNPATH_SUPPORT="no"
CONFIG_RUNPATH_FLAGS=""
dnl if test x"$ac_enable_new_dtags_ok" = x"yes"; then
if [[ "x${ac_enable_new_dtags_ok}" == "xyes" ]]; then
   HAS_RUNPATH_SUPPORT="yes"
   CONFIG_RUNPATH_FLAGS="-Wl,--enable-new-dtags"
fi
AC_SUBST(HAS_RUNPATH_SUPPORT)
AC_SUBST(CONFIG_RUNPATH_FLAGS)

dnl-------------------------------------------------------------------
dnl check if demangling using GNU C++ ABI is supported
dnl-------------------------------------------------------------------
dnl only do it if we have execinfo!
if [[ "x${ac_cv_header_execinfo_h}" == "xyes" ]]; then
    dnl check if demangling is enabled
    AC_ARG_ENABLE(demangling,
        [  --enable-demangling     enables demangling in backtrace (if supported) [default=no]],
        [enable_demangling_flag="$enableval"],
        [enable_demangling_flag="no"])

    dnl if asked for, check if deangling is supported
    if [[ "x${enable_demangling_flag}" == "xyes" ]] ; then
        ACX_CHECK_DEMANGLE_SUPPORT([AC_DEFINE(HAVE_DEMANGLING_SUPPORT, [], [defined if demangling is enabled at configure time and is supported through the GNU C++ ABI])],[])
    fi
fi

dnl substitutions
AC_SUBST(AM_CXXFLAGS)

dnl create the substitutions variables for the configuration script
AC_SUBST(PACKAGE_VERSION)
AC_SUBST(HAS_SHARED)

dnl-------------------------------------------------------------------
dnl plugin section
dnl
dnl each plugin can be enabled or disabled
dnl We also define a ENABLED_PLUGIN_<name> in config.h when plugin
dnl <name> is enabled
dnl WARNING: for make dist reasons all Makefiles MUST be created
dnl
dnl See m4/plugins.m4 to see what is enabled by default
dnl-------------------------------------------------------------------
dnl general
AC_CONFIG_FILES( plugins/Makefile )
CONFIG_LIBS_PLUGINS=""
CONFIG_LIBS_PLUGINS_STATIC=""

dnl deal with a couple of arguments to allow a whole bunch of plugins at a time
ACX_INIT_PLUGINS
ACX_CHECK_ALLPLUGINS
ACX_CHECK_ALLCXXPLUGINS

dnl initialise the list of allowed plugins
dnl This is meant to enable fastjet-config --list-plugins
LIST_ALL_PLUGINS=""

dnl allow to build all plugins in the same lib
AC_ARG_ENABLE(monolithic,
[  --enable-monolithic  Build all the (compiled) plugins in a single lib [default=yes]],
[case "${enableval}" in
  yes) monolithic=yes ;;
  no)  monolithic=no ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-monolithic) ;;
esac],[monolithic=yes])
AM_CONDITIONAL(MONOLITHIC_PLUGINS, test x$monolithic = xyes)
if [[ "x$monolithic" == "xyes" ]] ; then
   CONFIG_LIBS_PLUGINS=" -lfastjetplugins "
   CONFIG_LIBS_PLUGINS_STATIC=" \${installationdir}/lib/libfastjetplugins.a "
fi

dnl SISCone
ACX_CHECK_PLUGIN(SISCone, siscone, SISCONE, yes,
	         [CONFIG_LIBS_PLUGINS=${CONFIG_LIBS_PLUGINS}" -lsiscone_spherical -lsiscone"
                  CONFIG_LIBS_PLUGINS_STATIC=${CONFIG_LIBS_PLUGINS_STATIC}" \${installationdir}/lib/libsiscone.a \${installationdir}/lib/libsiscone_spherical.a"])
AC_CONFIG_SUBDIRS(plugins/SISCone/siscone)
dnl alternatively, to build SISCone directly from the FastJet build
dnl system, comment out the previoous lines and uncomment the ones
dnl below. Note that the main drawback of that method is that the
dnl SISCone version number need to be hardcoded in config.h.siscone.in
dnl   AM_CONFIG_HEADER(config.siscone.in.h:plugins/SISCone/siscone/siscone/config.h.in)
dnl   AC_CONFIG_FILES(plugins/SISCone/siscone/Makefile plugins/SISCone/siscone/siscone/Makefile plugins/SISCone/siscone/siscone/spherical/Makefile plugins/SISCone/siscone/examples/Makefile plugins/SISCone/siscone/examples/events/Makefile)

dnl CDFCones
ACX_CHECK_PLUGIN(CDFCones, cdfcones, CDFCONES, yes)
AC_CONFIG_FILES(plugins/CDFCones/CDFcode/Makefile)

dnl PxCone
ACX_CHECK_PLUGIN(PxCone, pxcone, PXCONE, no)

dnl D0RunIICone
ACX_CHECK_PLUGIN(D0RunIICone, d0runiicone, D0RUNIICONE, no)

dnl NestedDefs
ACX_CHECK_PLUGIN(NestedDefs, nesteddefs, NESTEDDEFS, yes)

dnl TrackJet
ACX_CHECK_PLUGIN(TrackJet, trackjet, TRACKJET, no)

dnl ATLASCone
ACX_CHECK_PLUGIN(ATLASCone, atlascone, ATLASCONE, no)

dnl CMSIterativeCone
ACX_CHECK_PLUGIN(CMSIterativeCone, cmsiterativecone, CMSITERATIVECONE, no)

dnl EECambridge
ACX_CHECK_PLUGIN(EECambridge, eecambridge, EECAMBRIDGE, yes)

dnl Jade
ACX_CHECK_PLUGIN(Jade, jade, JADE, yes)

dnl D0RunICone
ACX_CHECK_PLUGIN(D0RunICone, d0runicone, D0RUNICONE, no)

dnl GridJet
ACX_CHECK_PLUGIN(GridJet, gridjet, GRIDJET, yes)

dnl end of plugin section --- don't modify this line

dnl-------------------------------------------------------------------
dnl enable libtool (has to be done after the previous lib checks!)
dnl-------------------------------------------------------------------
AC_PROG_LIBTOOL

dnl-------------------------------------------------------------------
dnl add Fortran flags if any plugin that needs it is enabled
dnl-------------------------------------------------------------------
if [[ "x$ENABLE_PXCONE" == "xyes" ]] ; then
   dnl set up the fortran compiler and libs
   if test "x$FCFLAGS" == "x"; then
      FCFLAGS="-O2"
   fi
   ## Try to respect users' Fortran compiler variables
   if test "x$FC" == "x"; then
      if test "x$F77" != "x"; then
      	 FC="$F77"
      fi
   fi
   AC_PROG_FC
   if test "x$FC" == "x"; then
      AC_MSG_ERROR([Fortran compiler required to build PxCone but not found or not functional])
   fi
   AC_FC_LIBRARY_LDFLAGS
   F77=$FC
   AM_FFLAGS=$AM_FCFLAGS
   FFLAGS="$FFLAGS $FCFLAGS"
   FLIBS="$FLIBS $FCLIBS"
   AC_SUBST(FCFLAGS)
   AC_SUBST(FFLAGS)
   AC_SUBST(FLIBS)
   AC_SUBST(F77)
   AC_SUBST(AM_FFLAGS)
   AC_SUBST(AM_FCFLAGS)

   CONFIG_LIBS_PLUGINS=${CONFIG_LIBS_PLUGINS}" "${FLIBS}
   CONFIG_LIBS_PLUGINS_STATIC=${CONFIG_LIBS_PLUGINS_STATIC}" "${FLIBS}
fi

dnl-------------------------------------------------------------------
dnl be a bit clearer and use only AM_CXXFLAGS
dnl WARNING: this it totally unnecessary! Except for the D0 plugin:
dnl          since the CXXFLAGS are placed after the AM_CXXFLAGS in
dnl          the command line, a potential -O0 would be spoiled by
dnl          the -O3 in CXXFLAGS. By putting  CXXFLAGS in AM_CXXFLAGS
dnl          we have a better control on the order
dnl------------------------------------------------------------------
if [[ "x$USE_CUSTOM_FLAGS" == "xno" ]]; then
   AM_CXXFLAGS=${CXXFLAGS}" "${AM_CXXFLAGS}
   CXXFLAGS=""
fi


dnl------------------------------------------------------------------
dnl there is a problem with optimisation of the D0 code.
dnl One of the solutions to this issue is to disble inlie functions
dnl We ony do that if the compiler supports the option, otherwise, 
dnl we disable optimisation for that plugin
dnl------------------------------------------------------------------
if [[ "x$ENABLE_D0RUNIICONE" == "xyes" ]] ; then
   # check if the compiler supports -fno-inline
   AC_MSG_CHECKING([if $CXX supports -fno-inline])
   SAVE_CXXFLAGS="$CXXFLAGS"
   CXXFLAGS="-Werror -fno-inline"
   AC_LANG_PUSH(C++)
   AC_LINK_IFELSE([AC_LANG_SOURCE([int main(void){ return 0;}])],
      [ac_support_no_inline=yes],[ac_support_no_inline=no])
   AC_LANG_POP(C++)
   CXXFLAGS="$SAVE_CXXFLAGS"
   AC_MSG_RESULT([$ac_support_no_inline])
   D0_CXXFLAGS=""

   dnl if test x"$ac_enable_new_dtags_ok" = x"yes"; then
   if [[ "x${ac_support_no_inline}" == "xyes" ]]; then
      D0_CXXFLAGS=${AM_CXXFLAGS}" -fno-inline"
   else
      # in case of multiple -O options, the last one prevails
      D0_CXXFLAGS=${AM_CXXFLAGS}" -O0"
   fi
   
   AC_SUBST(D0_CXXFLAGS)
fi

dnl-------------------------------------------------------------------
dnl substitutions useful to know what we have to do in fastjet-config
dnl-------------------------------------------------------------------
AC_SUBST(CONFIG_LIBS_PLUGINS)
AC_SUBST(CONFIG_LIBS_PLUGINS_STATIC)
AC_SUBST(CONFIG_CXXFLAGS)
AC_SUBST(CONFIG_LIBS)
AC_SUBST(CONFIG_RUNPATH_FLAGS)
AC_SUBST(LIST_ALL_PLUGINS)

dnl-------------------------------------------------------------------
dnl tools section
dnl-------------------------------------------------------------------
AC_CONFIG_FILES( tools/Makefile tools/fastjet/Makefile tools/fastjet/tools/Makefile )


dnl-------------------------------------------------------------------
dnl make sure examples dir is there
dnl-------------------------------------------------------------------
AC_CONFIG_FILES( example/Makefile example/python/Makefile example/data/Makefile example/root/Makefile )

dnl-------------------------------------------------------------------
dnl and the documentation
dnl-------------------------------------------------------------------
AC_CONFIG_FILES( doc/Makefile doc/figs/Makefile )

dnl-------------------------------------------------------------------
dnl configuration script
dnl-------------------------------------------------------------------
AC_CONFIG_FILES([fastjet-config:fastjet-config.in], [chmod u+x fastjet-config])

dnl-------------------------------------------------------------------
dnl sample fortran wrapper
dnl-------------------------------------------------------------------
AC_CONFIG_FILES( fortran_wrapper/Makefile )

dnl-------------------------------------------------------------------
dnl python interface
dnl-------------------------------------------------------------------
AC_MSG_CHECKING([whether to build the python interface])
AC_ARG_ENABLE(pyext,
	[  --enable-pyext          Enable the python interface [default=no]],
	[ENABLE_PYEXT_FLAG="$enableval"],
	[ENABLE_PYEXT_FLAG="no"])
AC_MSG_RESULT(${ENABLE_PYEXT_FLAG})
if [[ "$ENABLE_PYEXT_FLAG" == "yes" ]] ; then
    AM_PATH_PYTHON
    AC_ARG_VAR([PYTHON_INCLUDE], [Include flags for python, bypassing python-config])
    AC_ARG_VAR([PYTHON_CONFIG], [Path to python-config])
    AS_IF([test -z "$PYTHON_INCLUDE"], [
        AS_IF([test -z "$PYTHON_CONFIG"], [
            AC_PATH_PROGS([PYTHON_CONFIG],
                          [python$PYTHON_VERSION-config python-config],
                          [no],
                          [`dirname $PYTHON`])
            AS_IF([test "$PYTHON_CONFIG" = no], [AC_MSG_ERROR([cannot find python-config for $PYTHON.])])
        ])
        AC_MSG_CHECKING([python include flags])
        PYTHON_INCLUDE=`$PYTHON_CONFIG --includes`
        AC_MSG_RESULT([$PYTHON_INCLUDE])
    ])

    # at a deeper level, one may want to be able to run swig dynamically
    #   
    # In the end we want to ship swig_wrapper.cpp and fastjet.py together
    # with the distribution but for development purposes (and possile
    # incoherences w end-user configurations) it is good to have teh
    # option to run swig automatically
    AC_MSG_CHECKING([whether to run swig])
    AC_ARG_ENABLE(swig,
                  [  --enable-swig          Runs swig to update the python interface [default=no]],
                  [ENABLE_SWIG_FLAG="$enableval"],
                  [ENABLE_SWIG_FLAG="no"])
    AC_MSG_RESULT(${ENABLE_SWIG_FLAG})
    if [[ "$ENABLE_SWIG_FLAG" == "yes" ]] ; then
        AX_PKG_SWIG([3.0.0], [], [AC_MSG_ERROR([SWIG explicitly requested but not found...]) ])
        SWIG_ENABLE_CXX
        SWIG_PYTHON
    fi
else
    ENABLE_SWIG_FLAG="no"
fi
AM_CONDITIONAL(BUILD_PYEXT, test x$ENABLE_PYEXT_FLAG = xyes)
AM_CONDITIONAL(BUILD_INCLUDES_SWIG, test x$ENABLE_SWIG_FLAG = xyes)
AC_CONFIG_FILES( pyinterface/Makefile )

dnl-------------------------------------------------------------------
dnl split the version number in major.minor.patch[prerelease] version
dnl
dnl We accordingly define FASTJET_VERSION_{MAJOR,MINOR,PATCHLEVEL} as
dnl well as FASTJET_VERSION_NUMBER. The latter is a 5-digit integer of
dnl the form 'XYYZZ' (for v X.Y.Z) that can be used to do version
dnl number comparisons.
dnl
dnl We assume that the 3 numbers are present + an optional pre-release
dnl tag. If the pre-release tag starts with a -, that - will be
dnl removed. All version numbers should only be made of numbers.
dnl
dnl Note that there are multiple strategies one can adopt for
dnl treatment of the pre-release tag (if any). Here we followed the
dnl idea that, being a candidate release, pre-releases version number
dnl will behave as their corresponiding stsable release. E.g. testing
dnl that the version is at least 3.1.0 using
dnl   #if FASTJET_VERSION_NUMBER >= 30100
dnl     ...
dnl   #endif
dnl will give a positive anwser for FastJet 3.1.0a1.
dnl-------------------------------------------------------------------
versionmajor=`echo $PACKAGE_VERSION | sed 's/^\([[^\.]]*\)\.\([[^\.]]*\)\.\([[^\.]]*\).*/\1/'`
versionminor=`echo $PACKAGE_VERSION | sed 's/^\([[^\.]]*\)\.\([[^\.]]*\)\.\([[^\.]]*\).*/\2/'`
versionending=`echo $PACKAGE_VERSION | sed 's/^\([[^\.]]*\)\.\([[^\.]]*\)\.\(.*\)$/\3/'`
versionpatchlevel=`echo $versionending | sed 's/^\([[0-9]][[0-9]]*\).*/\1/'`
versionprerelease=`echo $versionending | sed 's/^[[0-9-]]*//'`
versionnumber=`printf "%d%.2d%.2d" "$versionmajor" "$versionminor" "$versionpatchlevel"`

AC_DEFINE_UNQUOTED([VERSION_MAJOR],      [$versionmajor],      [Major version of this package])
AC_DEFINE_UNQUOTED([VERSION_MINOR],      [$versionminor],      [Minor version of this package])
AC_DEFINE_UNQUOTED([VERSION_PATCHLEVEL], [$versionpatchlevel], [Patch version of this package])
if [[ "x$versionprerelease" != "x" ]]; then
   AC_DEFINE_UNQUOTED([VERSION_PRERELEASE], ["$versionprerelease"], [Pre-release version of this package])
fi
AC_DEFINE_UNQUOTED([VERSION_NUMBER], [$versionnumber], [Version of the package under the form XYYZZ (instead of X.Y.Z)])

dnl AC_DEFINE_UNQUOTED([VERSION_IS(X,Y,Z)],
dnl  [(FASTJET_VERSION_MAJOR == (X) && FASTJET_VERSION_MINOR == (Y) && FASTJET_VERSION_PATCHLEVEL == (Z))],
dnl  [macro to require a precise version of FastJet (only works for FJ>=3.0.1)])
dnl 
dnl AC_DEFINE([VERSION_IS_AT_LEAST(X,Y,Z)],
dnl  [(FASTJET_VERSION_MAJOR > (X) || (FASTJET_VERSION_MAJOR == (X) && (FASTJET_VERSION_MINOR > (Y) || (FASTJET_VERSION_MINOR == (Y) && (FASTJET_VERSION_PATCHLEVEL >= (Z))))))],
dnl  [macro to require a minimal version of FastJet (only works for FJ>=3.0.1)])

dnl-------------------------------------------------------------------
dnl main library outputs
dnl-------------------------------------------------------------------
AC_CONFIG_FILES( Makefile src/Makefile include/Makefile include/fastjet/Makefile include/fastjet/internal/Makefile )

dnl-------------------------------------------------------------------
dnl give a nice summary of the config
dnl-------------------------------------------------------------------
dnl echo ""
dnl echo "Configuration summary:"
dnl echo "----------------------"
dnl echo "  Installation directory     "${prefix}
dnl echo "  Shared libraries           "${enable_shared}
dnl echo "  Static libraries           "${enable_static}
dnl echo "  Debug flag                 "$ENABLE_DEBUG_FLAG
dnl echo "  CGAL support               "${ENABLE_CGAL_FLAG}
dnl AC_PLUGIN_SUMMARY
dnl echo "  Monolithic plugins lib     "$monolithic

CONFIG_SUMMARY=""
CONFIG_SUMMARY+="Configuration summary:\n"
CONFIG_SUMMARY+="----------------------\n"
if [[ "x${prefix}" == "xNONE" ]] ; then
CONFIG_SUMMARY+="  Installation directory     /usr/local\n"
else
CONFIG_SUMMARY+="  Installation directory     "${prefix}"\n"
fi
CONFIG_SUMMARY+="  Shared libraries           "${enable_shared}"\n"
CONFIG_SUMMARY+="  Static libraries           "${enable_static}"\n"
CONFIG_SUMMARY+="  Debug flag                 "$ENABLE_DEBUG_FLAG"\n"
CONFIG_SUMMARY+="  CGAL support               "${ENABLE_CGAL_FLAG}"\n"
if [[ "$ENABLE_THREAD_SAFETY_FLAG" == "yes" ]]; then
CONFIG_SUMMARY+="  Thread safety              yes (requires C++11)\n"
else
if [[ "$ENABLE_LIMITED_THREAD_SAFETY_FLAG" == "yes" ]]; then
CONFIG_SUMMARY+="  Thread safety              limited (requires C++11)\n"
else
CONFIG_SUMMARY+="  Thread safety              no\n"
fi
fi
AC_PLUGIN_SUMMARY([yes])
CONFIG_SUMMARY+="  Monolithic plugins lib     "$monolithic"\n"
CONFIG_SUMMARY+="  Python interface           "${ENABLE_PYEXT_FLAG}"\n"
CONFIG_SUMMARY+="  Swig python generator      "${ENABLE_SWIG_FLAG}"\n"
AC_SUBST(CONFIG_SUMMARY)


CONFIGURE_INVOCATION="$0 $ac_configure_args"
AC_SUBST(CONFIGURE_INVOCATION)

dnl-------------------------------------------------------------------
dnl finish with a nice message
dnl-------------------------------------------------------------------
AC_OUTPUT
printf "$CONFIG_SUMMARY"
echo "---------------------------------------------------"
echo "        IMPORTANT NOTICE TO FASTJET USERS          "
echo "                                                   "
echo " When using the FastJet library, you can retrieve  "
echo " compilation and link flags using fastjet-config.  "
echo " Type 'fastjet-config --help' for more information "
echo "---------------------------------------------------"

