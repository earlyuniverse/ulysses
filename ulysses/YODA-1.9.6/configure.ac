## Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([YODA],[1.9.6],[yoda@projects.hepforge.org],[YODA])

## Check and block installation into the src/build dir
if test "$prefix" = "$PWD"; then
  AC_MSG_ERROR([Installation into the build directory is not supported: use a different --prefix argument])
fi
## Force default prefix to have a path value rather than NONE
if test "$prefix" = "NONE"; then
   prefix=/usr/local
fi

AC_CONFIG_SRCDIR([src/Counter.cc])
AM_INIT_AUTOMAKE([subdir-objects -Wall dist-bzip2 1.10])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_HEADERS([include/YODA/Config/DummyConfig.h include/YODA/Config/YodaConfig.h include/YODA/Config/BuildConfig.h])
AC_DEFINE_UNQUOTED(YODA_VERSION, "$PACKAGE_VERSION", "YODA version string")
AC_DEFINE_UNQUOTED(YODA_NAME, "$PACKAGE_NAME", "YODA name string")
AC_DEFINE_UNQUOTED(YODA_STRING, "$PACKAGE_STRING", "YODA name and version string")
AC_DEFINE_UNQUOTED(YODA_TARNAME, "$PACKAGE_TARNAME", "YODA short name string")
AC_DEFINE_UNQUOTED(YODA_BUGREPORT, "$PACKAGE_BUGREPORT", "YODA contact email address")

## OS X
AC_CEDAR_OSX

## Set default compiler flags
if test "x$CXXFLAGS" = "x"; then CXXFLAGS="-O3"; fi

## Compiler setup
AC_LANG(C++)
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX([11], [noext], [mandatory])

AC_PROG_INSTALL
AC_PROG_LN_S
AC_DISABLE_STATIC
AC_PROG_LIBTOOL

## Work out library suffix for the build
LIB_SUFFIX=\\\"$shrext_cmds\\\"
AC_SUBST([LIB_SUFFIX])

## Set default build flags
AC_CEDAR_CHECKCXXFLAG([-pedantic], [AM_CXXFLAGS="$AM_CXXFLAGS -pedantic"])
AC_CEDAR_CHECKCXXFLAG([-Wall], [AM_CXXFLAGS="$AM_CXXFLAGS -Wall -Wno-format"])
dnl AC_CEDAR_CHECKCXXFLAG([-std=c++98], [AM_CXXFLAGS="$AM_CXXFLAGS -std=c++98"])
dnl AC_CEDAR_CHECKCXXFLAG([-Wno-unused-variable], [AM_CXXFLAGS="$AM_CXXFLAGS -Wno-unused-variable"])

## Debug flag (default=none)
AC_ARG_ENABLE([debug], [AC_HELP_STRING(--enable-debug, [build with debugging symbols  @<:@default=no@:>@])], [], [enable_debug=no])
if test x$enable_debug = xyes; then
  [AM_CXXFLAGS="$AM_CXXFLAGS -g"]
fi


## Optional zlib support for gzip-compressed data streams/files
AX_CHECK_ZLIB

## Optional ROOT compatibility
AC_ARG_ENABLE([root], [AC_HELP_STRING(--disable-root,
  [don't try to build YODA interface to PyROOT (needs root-config) @<:@default=yes@:>@])], [], [enable_root=yes])
if test "x$enable_root" = "xyes"; then
  AC_PATH_PROG(ROOTCONFIG, [root-config])
  if test "x$ROOTCONFIG" = "x"; then
    enable_root=no;
    AC_MSG_WARN([root-config not found -- not building extra ROOT compatibility tools])
  else
    AC_MSG_CHECKING([ROOT version])
    ROOT_VERSION=`$ROOTCONFIG --version`
    ROOT_MAJOR_VERSION=`echo $ROOT_VERSION | cut -d. -f1`
    ROOT_MINOR_VERSION=`echo $ROOT_VERSION | cut -d. -f2 | cut -d/ -f1`
    ROOT_MICRO_VERSION=`echo $ROOT_VERSION | cut -d. -f2 | cut -d/ -f2`
    AC_MSG_RESULT([$ROOT_VERSION ($ROOT_MAJOR_VERSION,$ROOT_MINOR_VERSION,$ROOT_MICRO_VERSION)])
    if test "$ROOT_MAJOR_VERSION" -lt 6; then
      enable_root=no;
      AC_MSG_WARN([ROOT major version is < 6 -- not building extra ROOT compatibility tools])
    fi
    # TODO: Test for existence of TPython, instance_from_void API, etc.
    #AM_CXXFLAGS="$AM_CXXFLAGS -Wno-long-long"
    ROOT_CXXFLAGS=`$ROOTCONFIG --cflags`
    ROOT_LDFLAGS=`$ROOTCONFIG --ldflags`
    ROOT_LIBS=`$ROOTCONFIG --libs`
    if test "$ROOT_MAJOR_VERSION" -eq 6 && test "$ROOT_MINOR_VERSION" -lt 22; then
      ROOT_LIBS="$ROOT_LIBS -lPyROOT"
    else
      ROOT_LIBS="$ROOT_LIBS -lROOTTPython"
    fi
    AC_SUBST(ROOT_CXXFLAGS)
    AC_SUBST(ROOT_LDFLAGS)
    AC_SUBST(ROOT_LIBS)
  fi
fi
AM_CONDITIONAL(ENABLE_ROOT, [test x$enable_root = xyes])
if test x$enable_root = xyes; then
  AC_MSG_NOTICE([Building extra ROOT compatibility tools])
else
  AC_MSG_NOTICE([Not building extra ROOT compatibility tools])
fi


## Python extension
AC_ARG_ENABLE(pyext, [AC_HELP_STRING(--disable-pyext,
  [don't build Python module (default=build)])],
  [], [enable_pyext=yes])
## Basic Python checks
if test x$enable_pyext = xyes; then
  AX_PYTHON_DEVEL([>= '2.7.3'])
  PYTHON_FULL_VERSION=`$PYTHON -c 'from __future__ import print_function; import platform; print(platform.python_version())'`
  PYTHON_MAJOR_VERSION=`$PYTHON -c 'from __future__ import print_function; import sys; print(sys.version_info.major)'`
  AC_SUBST(PYTHON_VERSION)
  AC_SUBST(PYTHON_FULL_VERSION)
  AC_SUBST(PYTHON_MAJOR_VERSION)
  YODA_PYTHONPATH=`$PYTHON -c "from __future__ import print_function; import distutils.sysconfig; print(distutils.sysconfig.get_python_lib(prefix='$prefix', plat_specific=True));"`
  AC_SUBST(YODA_PYTHONPATH)
  if test -z "$PYTHON"; then
    AC_MSG_ERROR([Can't build Python extension since python can't be found])
    enable_pyext=no
  fi
  if test -z "$PYTHON_CPPFLAGS"; then
    AC_MSG_ERROR([Can't build Python extension since Python.h header file cannot be found])
    enable_pyext=no
  fi
fi
AM_CONDITIONAL(ENABLE_PYEXT, [test x$enable_pyext == xyes])


## Cython checks
if test x$enable_pyext == xyes; then
  AM_CHECK_CYTHON([0.24], [:], [:])
  if test x$CYTHON_FOUND = xyes; then
    AC_PATH_PROGS(CYTHON, [$CYTHON cython-$PYTHON_VERSION cython$PYTHON_VERSION cython-$PYTHON_MAJOR_VERSION cython$PYTHON_MAJOR_VERSION cython])
    if test x$CYTHON != x; then
        AC_MSG_NOTICE([Cython >= 0.24 found at $CYTHON: Python extension source can be rebuilt (for developers)])
        # Force rebuild since we have a sufficient Cython
        for i in pyext/yoda//*.pyx; do touch $i; done
    fi
  fi

  AC_CHECK_FILE([pyext/yoda/core.cpp],
                [],
                [if test "x$CYTHON_FOUND" != "xyes"; then
                  AC_MSG_ERROR([Cython is required for --enable-pyext, no pre-built core.cpp was found.])
                fi])

  ## Set extra Python extension build flags (to cope with Cython output code oddities)
  PYEXT_CXXFLAGS=$CXXFLAGS
  AC_CEDAR_CHECKCXXFLAG([-Wno-unused-but-set-variable], [PYEXT_CXXFLAGS="$PYEXT_CXXFLAGS -Wno-unused-but-set-variable"])
  AC_CEDAR_CHECKCXXFLAG([-Wno-sign-compare], [PYEXT_CXXFLAGS="$PYEXT_CXXFLAGS -Wno-sign-compare"])
  AC_CEDAR_CHECKCXXFLAG([-Wno-strict-prototypes], [PYEXT_CXXFLAGS="$PYEXT_CXXFLAGS -Wno-strict-prototypes"])
  AC_SUBST(PYEXT_CXXFLAGS)
  AC_MSG_NOTICE([All Python build checks successful: 'yoda' Python extension will be built])
fi
AM_CONDITIONAL(WITH_CYTHON, [test x$CYTHON_FOUND = xyes])


## Extend and substitute the default build flags after lib testing
AM_CPPFLAGS="-I\$(top_srcdir)/include -I\$(top_builddir)/include"
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_CXXFLAGS)

## setup.py puts its build artifacts into a labelled path
## this helps the test scripts to find them locally instead of
## having to install first
YODA_SETUP_PY_PATH=$(${PYTHON} -c 'from __future__ import print_function; import distutils.util as u, sys; vi=sys.version_info; print("lib.%s-%s.%s" % (u.get_platform(),vi.major, vi.minor))')
AC_SUBST(YODA_SETUP_PY_PATH)

## Build Doxygen if possible
AC_PATH_PROG(DOXYGEN, doxygen)
AM_CONDITIONAL(WITH_DOXYGEN, test "$DOXYGEN")


## Build file output
AC_EMPTY_SUBST
AC_CONFIG_FILES([Makefile Doxyfile])
AC_CONFIG_FILES([include/Makefile include/YODA/Makefile])
AC_CONFIG_FILES([src/Makefile
                 src/tinyxml/Makefile
                 src/yamlcpp/Makefile
])
AC_CONFIG_FILES([tests/Makefile])
AC_CONFIG_FILES([pyext/Makefile
                 pyext/setup.py
                 pyext/yoda/Makefile
                 ])
AC_CONFIG_FILES([bin/Makefile bin/yoda-config])
AC_CONFIG_FILES([yodaenv.sh yoda.pc])

AC_OUTPUT

if test x$enable_pyext == xyes; then
   cat <<EOF

************************************************************
YODA CONFIGURED!

Now build and install (to the $prefix tree) with e.g.
make -j2 && make -j2 install
************************************************************
EOF
fi
