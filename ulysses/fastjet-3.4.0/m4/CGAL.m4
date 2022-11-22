dnl CHECK CGAL BEGIN
dnl This script takes too arguments, the action on success and the action on failure.
dnl It first checks if a --with-cgalmakefile) value is passed.
dnl If that fails, it checks if the CGAL_MAKEFILE environment variable
dnl is defined. If that fails too, it seaches for CGAL in the standard places.
dnl
dnl If a Makefile is specified and not found, an error is returned
dnl If the standard Makefile is not found, then one check for the 
dnl   --with-cgaldir presence and check if CGAL can be found in standard places
dnl 
dnl Note that in the last case, we use a predefined set of compiler and linker
dnl arguments, a method that may fail.
dnl
dnl CGAL_CXXFLAGS, CGAL_CPPFLAGS, CGAL_LDFLAGS and CGAL_LIBS are all defined.
AC_DEFUN([ACX_CHECK_CGAL],
[
dnl Start by looking if a Makefile is given through --with-cgalmakefile
acx_cgal_found=no
acx_cgal_default_makefile=no

AC_ARG_WITH(cgalmakefile,
            [AC_HELP_STRING([--with-cgalmakefile=makefile], [Use the following CGAL makefile (CGAL <= 3.3.x)])])

dnl ckeck if a directory is specified for CGAL
dnl Note: this is of no effect if a Makefile is used
AC_ARG_WITH(cgaldir,
            [AC_HELP_STRING([--with-cgaldir=dir], [Assume the given directory for CGAL (CGAL >= 3.4)])])

dnl allows to specify a non-standard installation directory for Boost used by CGAL
AC_ARG_WITH(cgal_boostdir,
            [AC_HELP_STRING([--with-cgal-boostdir=dir], [Assume the given directory for Boost needed by CGAL (CGAL >= 3.4; requires an installed Boost)])])

dnl allows to specify a non-standard installation directory for GMP used by CGAL
AC_ARG_WITH(cgal_gmpdir,
            [AC_HELP_STRING([--with-cgal-gmpdir=dir], [Assume the given directory for GMP needed by CGAL])])

dnl allows to specify a non-standard installation directory for MPFR used by CGAL
AC_ARG_WITH(cgal_mpfrdir,
            [AC_HELP_STRING([--with-cgal-mpfrdir=dir], [Assume the given directory for MPFR needed by CGAL])])

dnl allows one to say that the CGAL build is header-only
AC_ARG_ENABLE(cgal-header-only,
              [  --enable-cgal-header-only   enable build with header-only install of CGAL, e.g. as for CGALv5; in that case do not use --enable-cgal [default=no]],
              [ENABLE_CGAL_HEADER_ONLY_FLAG="$enableval"],
              [ENABLE_CGAL_HEADER_ONLY_FLAG="no"])

dnl ------------------------------------------------------------------------
dnl define CGAL_MAKEFILE to be 
dnl  1. the value given to --with-cgalmakefile
dnl  2. the environment var
dnl  3. the default location: /usr/share/cgal/cgal.mk
dnl in that order 
dnl Note: using --with-cgalmakefile=no disables CGAL. In our case this is redundant as 
dnl       we already have a --enable-cgal option, but we keep it as it was in the
dnl       original code
case $with_cgalmakefile in
    yes | "")
        dnl check if a env. var is defined
        if test x${CGAL_MAKEFILE} = x ; then
            # if no installation dir is explicitly specified, try a default makefile
            if test \! -z "$with_cgaldir"; then
                CGAL_MAKEFILE=""
            else 
                CGAL_MAKEFILE="/usr/share/cgal/cgal.mk" 
            fi
        fi
        acx_cgal_default_makefile=yes;
        ;;
    no) acx_cgal_found=disable ;;
    -* | */* | *.a | *.so | *.so.* | *.o) CGAL_MAKEFILE="$with_cgalmakefile" ;;
    *) CGAL_MAKEFILE="$with_cgalmakefile" ;;
esac

dnl ------------------------------------------------------------------------
dnl check for the presence of the CGAL Makefile
dnl if present, get the compilation and linker flags there
if test "$acx_cgal_found" == no; then
    AC_MSG_CHECKING(CGAL_MAKEFILE)

    dnl check that a Makefile has been specified
    if test \! -z "$CGAL_MAKEFILE"; then

        dnl make sure the file exists
        if test -e "$CGAL_MAKEFILE"; then
            dnl to get the arguments, we'll just write a temporary Makefile
            dnl that "echo"s the CGAL_CXXFLAGS and CGAL_LDFLAGS variables upon request
            tname=`mktemp /tmp/cgal_makefile_dsrXXXXXX`
    
cat > $tname << _ACEOF
include $CGAL_MAKEFILE

cppflags:
	@echo \$(CGAL_CXXFLAGS)

cxxflags:
	@echo
ldflags:
	@echo \$(CGAL_LDFLAGS)
_ACEOF

            dnl get the compiler and preprocessor info
            CGAL_CPPFLAGS="`make -s -f $tname cppflags`"
            CGAL_CXXFLAGS="`make -s -f $tname cxxflags`"
            
            dnl for the libs, put the -l directives in CGAL_LIBS and the others in 
            dnl CGAL_LDFLAGS
            CGAL_LDFLAGST="`make -s -f $tname ldflags`"
            for i in $CGAL_LDFLAGST; do
                if test `echo $i| grep -c ^-l`; then 
                    CGAL_LIBS="$CGAL_LIBS $i"
                else
                    CGAL_LDFLAGS="$CGAL_LDFLAGS $i"
                fi
            done

            dnl remove the temporary Makefile and tell we succeeded
            rm -f $tname
            AC_MSG_RESULT(yes)
            acx_cgal_found=yes
        else
            dnl the Makefile does not exists
            dnl We only need to return an error if the Makefile was not the default one
            AC_MSG_RESULT(inexistent)
            if test "$acx_cgal_default_makefile" == no; then
                AC_ERROR([CGAL_MAKEFILE defined as $CGAL_MAKEFILE, but the makefile does not exist.])
            fi
        fi
    else 
        dnl the Makefile is invalid
        AC_MSG_RESULT(not specified)
        dnl dnl We only need to return an error if the Makefile was not the default one
        dnl if test "$acx_cgal_default_makefile" == no; then
        dnl     AC_ERROR([CGAL_MAKEFILE defined as $CGAL_MAKEFILE, but the makefile does not exist.])
        dnl fi
    fi
fi

dnl ------------------------------------------------------------------------
dnl if CGAL hasn't been found yet and the Makefile search didn't return an error message
dnl we check if we can get the arguments ourselves. This method leave the option to specify
dnl a CGAL directory using the --with-cgaldir directive
dnl If no dir are specified, CGAL will be searched for in standard places
if test "$acx_cgal_found" == no; then

    dnl store the compiler/linker flags so we can restore them in case of failure
    dnl After each check we add found flags to both CGAL_... flags and global flags
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    save_CXXFLAGS="$CXXFLAGS"
    save_CPPFLAGS="$CPPFLAGS"

    dnl check support for the floating-point specifications needed for CGAL
    dnl  . gcc required -frounding-math
    dnl  . icc requires -fp-model srticy
    dnl    
    dnl Note that clang does not support it but simply issues a
    dnl warning. To avoid having that, we force it to be an error if it is
    dnl not supported
    dnl   
    dnl also, icpc supports -Werror but does not seem to convert
    dnl "unrecognised option" into an error, at least for icc 13.1.3
    dnl that I've tested. Let's live with that for hte time being since
    dnl it will just issue a $1006 command-line warning 
    ADDITIONAL_CGAL_FLAGS=""
    AC_LANG_PUSH(C++)
    AX_CHECK_COMPILER_FLAGS([-frounding-math -Werror],[ADDITIONAL_CGAL_FLAGS=${ADDITIONAL_CGAL_FLAGS}" -frounding-math"])
    AX_CHECK_COMPILER_FLAGS([-fp-model strict -Werror],[ADDITIONAL_CGAL_FLAGS=${ADDITIONAL_CGAL_FLAGS}" -fp-model strict"])
    AC_LANG_POP(C++)

    CGAL_CPPFLAGS="${ADDITIONAL_CGAL_FLAGS}"
    CGAL_LIBS=""

    dnl check if some path for CGAL is provided
    dnl
    dnl note that some platforms install using a "lib64" prefix, so
    dnl we'll use a stand-in that we'll replace later

    dnl First check if a CGAL installation dir has been supplied.
    if test \! -z "$with_cgaldir"; then
        AC_MSG_CHECKING(CGAL in ${with_cgaldir})
        if test \! -d ${with_cgaldir}/include; then
            AC_MSG_RESULT([no CGAL header directory found])
            CGAL_CPPFLAGS=""
            CGAL_LIBS=""
            $2
            exit
        fi
        CGAL_CPPFLAGS="$CGAL_CPPFLAGS -I${with_cgaldir}/include"
        libdir_found="no"
        for libdircandidate in lib lib64 lib32; do
            if test -d ${with_cgaldir}/${libdircandidate}; then
                libdir_found="yes"
                CGAL_LIBS="${CGAL_LIBS} -L${with_cgaldir}/${libdircandidate} -Wl,-rpath,${with_cgaldir}/${libdircandidate}"
                break
            fi
        done
        if test "$libdir_found" == yes; then
            AC_MSG_RESULT([${libdircandidate}])
        else 
            AC_MSG_RESULT([no CGAL library directory found])
            CGAL_CPPFLAGS=""
            CGAL_LIBS=""
            $2
            exit
        fi
    fi

    dnl if a non-standard Boost location has been specified, add it to
    dnl the compilation flags
    dnl
    dnl it seems that an explicit -lbbost is not needed to compile FJ
    dnl or codes using FJ. The code below just adds the -L flags in
    dnl case it helps as well as the include files
    dnl
    if test \! -z "$with_cgal_boostdir"; then
        AC_MSG_CHECKING(with Boost in ${with_cgal_boostdir})
        if test \! -d ${with_cgal_boostdir}/include; then
            AC_MSG_RESULT([no CGAL boost header directory found])
            CGAL_CPPFLAGS=""
            CGAL_LIBS=""
            $2
            exit
        fi
        CGAL_CPPFLAGS="$CGAL_CPPFLAGS -I${with_cgal_boostdir}/include"
        AC_MSG_RESULT([yes])
        dnl skip lib path search since it looks like FJ 
        dnl does not require the Boost libs to be functional
    fi
    dnl the following lines are commented out because it looks like FJ 
    dnl does not require the Boost libs to be functional
    dnl
    dnl AC_LANG_PUSH(C++)
    dnl AC_CHECK_HEADERS([boost/version.hpp])
    dnl AC_LANG_POP(C++)
    dnl AC_CHECK_LIB(boost_system, main, [CGAL_LIBS="-lboost_system $CGAL_LIBS"])
    dnl AC_CHECK_LIB(boost_thread, main, [CGAL_LIBS="-lboost_thread $CGAL_LIBS"])

    dnl check for gmp as CGAL depends on it
    dnl
    dnl it seems that FJ (and code using FJ) depend both on headers and lib
    dnl
    dnl if a non-standard location has been specified, add it to
    dnl the compilation flags
    if test \! -z "$with_cgal_gmpdir"; then
        AC_MSG_CHECKING(with GMP in ${with_cgal_gmpdir})
        if test \! -d ${with_cgal_gmpdir}/include; then
            AC_MSG_RESULT([no CGAL GMP header directory found])
            CGAL_CPPFLAGS=""
            CGAL_LIBS=""
            $2
            exit
        fi
        CGAL_CPPFLAGS="$CGAL_CPPFLAGS -I${with_cgal_gmpdir}/include"
        libdir_found="no"
        for libdircandidate in lib lib64 lib32; do
            if test -d ${with_cgal_gmpdir}/${libdircandidate}; then
                libdir_found="yes"
                CGAL_LIBS="${CGAL_LIBS} -L${with_cgal_gmpdir}/${libdircandidate} -Wl,-rpath,${with_cgal_gmpdir}/${libdircandidate}"
                break
            fi
        done
        if test "$libdir_found" == yes; then
            AC_MSG_RESULT([${libdircandidate}])
        else 
            AC_MSG_RESULT([no GMP lib found in directory])
            CGAL_CPPFLAGS=""
            CGAL_LIBS=""
            $2
            exit
        fi
    fi
    dnl search for gmp 
    CXXFLAGS=${save_CXXFLAGS}" $CGAL_CPPFLAGS"
    CPPFLAGS=${save_CPPFLAGS}" $CGAL_CPPFLAGS"
    LIBS="${save_LIBS}"
    LDFLAGS="$CGAL_LIBS ${save_LDFLAGS}"

    AC_CHECK_HEADERS(gmp.h)
    AC_CHECK_LIB(gmp, main, cgal_have_gmplib=yes, cgal_have_gmplib=no)
    dnl AC_CHECK_LIB(gmpxx, main, [CGAL_LIBS="$CGAL_LIBS -lgmpxx"])
    if test "$cgal_have_gmplib" == yes; then
        CGAL_LIBS="-lgmp $CGAL_LIBS"
    else
        CGAL_CPPFLAGS=""
        CGAL_LIBS=""
        LIBS="$save_LIBS"
        LDFLAGS="$save_LDFLAGS"
        CXXFLAGS="$save_CXXFLAGS"
        CPPFLAGS="$save_CPPFLAGS"
        $2
        exit
    fi        

    dnl check for MPFR as CGAL depends on it
    dnl
    dnl it seems that FJ (and code using FJ) depend only on the headers.
    dnl We include both on headers and lib
    dnl
    dnl if a non-standard location has been specified, add it to
    dnl the compilation flags
    dnl 
    dnl it seems that in the case of MPFR, only heades are needed
    if test \! -z "$with_cgal_mpfrdir"; then
        AC_MSG_CHECKING(with MPFR in ${with_cgal_mpfrdir})
        if test \! -d ${with_cgal_mpfrdir}/include; then
            AC_MSG_RESULT([no CGAL GMP header directory found])
            LIBS="$save_LIBS"
            LDFLAGS="$save_LDFLAGS"
            CXXFLAGS="$save_CXXFLAGS"
            CPPFLAGS="$save_CPPFLAGS"
            CGAL_CPPFLAGS=""
            CGAL_LIBS=""
            $2
            exit
        fi
        CGAL_CPPFLAGS="$CGAL_CPPFLAGS -I${with_cgal_mpfrdir}/include"
        AC_MSG_RESULT([yes])
    fi
    dnl search for mpfr
    CXXFLAGS=${save_CXXFLAGS}" $CGAL_CPPFLAGS"
    CPPFLAGS=${save_CPPFLAGS}" $CGAL_CPPFLAGS"
    LIBS="${save_LIBS}"
    LDFLAGS="$CGAL_LIBS ${save_LDFLAGS}"

    AC_CHECK_HEADERS(mpfr.h)
    dnl AC_CHECK_LIB(mpfr, main, [CGAL_LIBS="-lmpfr $CGAL_LIBS"])


    dnl search the CGAL headers
    dnl
    dnl set temporarily the check language to C++
    dnl The default choice otherwise in C whih is usually fine apart from 
    dnl the tests related to CGAL which have c-like header names and fail
    dnl with gcc. The switch to g++ allows CGAL to be detected if it is
    dnl not detected through a Makefile.
    dnl
    dnl Note that one can temporarily change lang. By using
    dnl   AC_LANG_PUSH(language)
    dnl followed by
    dnl   AC_LANG_POP(language)

    dnl first save the current flags so we can restore them later if needed

    CXXFLAGS=${save_CXXFLAGS}" $CGAL_CPPFLAGS"
    CPPFLAGS=${save_CPPFLAGS}" $CGAL_CPPFLAGS"

    AC_LANG_PUSH(C++)
    AC_CHECK_HEADER([CGAL/Exact_predicates_inexact_constructions_kernel.h], [cgal_have_header=yes], [cgal_have_header=no])
    AC_LANG_POP(C++)

    dnl if the headers have been found, check for the libs
    if test "$cgal_have_header" == yes; then
        if [[ "$ENABLE_CGAL_HEADER_ONLY_FLAG" == "no" ]] ; then
            dnl prepare the linker flags for test
            LIBS="${save_LIBS}"
            LDFLAGS="$CGAL_LIBS ${save_LDFLAGS}"
            
            AC_LANG_PUSH(C++)
            AC_CHECK_LIB(CGAL, main, cgal_have_lib=yes, cgal_have_lib=no)
            dnl AC_CHECK_LIB(CGALcore++, main, [CGAL_LIBS="$CGAL_LIBS -lCGALcore++"])
            AC_LANG_POP(C++)
            if test "$cgal_have_lib" == yes; then
                CGAL_LIBS=" -lCGAL "${CGAL_LIBS}
                
                dnl we can finally claim we've found CGAL!
                acx_cgal_found=yes
            fi
        else
            acx_cgal_found=yes
        fi
    fi 

    dnl restore the compiler flags to their default
    LIBS="$save_LIBS"
    LDFLAGS="$save_LDFLAGS"
    CXXFLAGS="$save_CXXFLAGS"
    CPPFLAGS="$save_CPPFLAGS"

    dnl if the lib has not been found, reset the saved vars to their original values
    if test "$acx_cgal_found" == no; then 
        CGAL_CPPFLAGS=""
        CGAL_LIBS=""
    fi
fi

dnl ------------------------------------------------------------------------
dnl sonclude whether CGAL has been found or not
AC_MSG_CHECKING(CGAL)
if test "$acx_cgal_found" == yes; then
    AC_MSG_RESULT(yes);
    $1
else
    AC_MSG_RESULT(no);
    $2
fi])

dnl CHECK CGAL END
