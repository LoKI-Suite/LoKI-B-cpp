AC_PREREQ([2.69])
AC_INIT([loki-b],
	[m4_esyscmd_s([cat VERSION.txt])],
	[loki@tecnico.ulisboa.pt])
AC_CONFIG_SRCDIR([Makefile.am])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([autodir])
AC_CONFIG_MACRO_DIR([autodir])
AM_INIT_AUTOMAKE([1.11 foreign subdir-objects])
AC_CONFIG_FILES([
	Makefile
	app/Makefile
	doc/Makefile
	doc/doxygen/Makefile
	doc/doxygen/loki-b.doxy
	lib/Makefile
	ideas/Makefile
	tests/Makefile
	])

LT_INIT([disable-static win32-dll])
AC_SUBST(LIBTOOL_DEPS)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

CXXFLAGS="$CXXFLAGS -g -pedantic -Wall"

dnl ***************************************************************************
dnl
dnl     C++ standard version
dnl
dnl Check that the compiler, with the provided or calculated compiler options,
dnl is running in C++2017 (or later) mode.
dnl
dnl ***************************************************************************

AC_LANG_PUSH(C++)

AC_TRY_COMPILE(
  [ ],
  [ static_assert(__cplusplus>=201703L, "No C++17 support."); ],
  [ AC_MSG_NOTICE([C++17 detected. No futher action needed.]) ],
  [
    dnl The C++17 test failed. Add the -stc++17 flag and try again.
    AC_MSG_WARN([No C++17 support. I will try again with compiler flag -std=c++17.])
    CXXFLAGS="$CXXFLAGS -std=c++17"
    AC_TRY_COMPILE(
      [ ],
      [ static_assert(__cplusplus>=201703L, "No C++17 support."); ],
      [ ],
      [
        AC_MSG_ERROR([C++17 or later is required. Reconfigure with appropriate CXXFLAGS, or update your compiler.]) ]
    )
  ]
)

dnl ***************************************************************************
dnl
dnl     The Eigen headers and libraries
dnl
dnl ***************************************************************************
AC_ARG_WITH(eigen,
	[  --with-eigen(=location of the eigen3 directory)
			  use the Eigen library. The default path is /usr/include/eigen3],
	with_eigen="${withval}", with_eigen="yes")
AC_MSG_CHECKING(--with-eigen)
AC_MSG_RESULT($with_eigen)
EIGEN_MATVEC_LIBS=
EIGEN_CPPFLAGS=
if test "$with_eigen" != "no" ; then
  dnl Get the eigen3 dir and add that to the include path
  if test "$with_eigen" = "yes" ; then
    dnl Assume eigen3 is in /usr/include/eigen3
    EIGEN_CPPFLAGS="-I/usr/include/eigen3"
  else
    EIGEN_CPPFLAGS="-I$with_eigen"
  fi
  CPPFLAGS="$CPPFLAGS $EIGEN_CPPFLAGS"
  dnl Check that we can find the Eigen headers in the specified location
  AC_CHECK_HEADER(Eigen/Core,, AC_MSG_ERROR(Eigen3 was not found in include path $EIGEN_CPPFLAGS.))
  dnl Subdirs for eigen-specific stuff
  EIGEN_SUBDIR="eigen"
  EIGEN_MATVEC_LIBS="$(pwd)/plmath/lib/eigen/libpleigen_lapack.la $LAPACK_LIBS"
  AC_DEFINE_UNQUOTED([HAVE_EIGEN],1,[Define if we have the Eigen matrix vector library])
fi

AC_SUBST(EIGEN_CPPFLAGS)

# nlohmann json
AC_ARG_WITH(json,
  [  --with-json(=json prefix)
              use the nlohmann json library],
  with_json="${withval}", with_json="yes")
AC_MSG_CHECKING(--with-json)
AC_MSG_RESULT($with_json)
JSON_CPPFLAGS=
if test "$with_json" != "no" ; then
  dnl Testsuites that rely on json
  if test "$with_json" = "yes" ; then
    AC_CHECK_HEADER(nlohmann/json.hpp,
      with_json=yes,
      AC_MSG_ERROR(The nlohann json headers were not found. Please install these and reconfigure.))
  else
    JSON_CPPFLAGS="-I$with_json"
    CPPFLAGS="$CPPFLAGS $JSON_CPPFLAGS"
    AC_CHECK_HEADER(nlohmann/json.hpp,
      with_json=yes,
      AC_MSG_ERROR(The nlohann json headers were not found. Please install these and reconfigure.))
  fi
  AC_DEFINE_UNQUOTED([HAVE_JSON],1,[Define if we have the nlohmann json library])
else
  AC_MSG_ERROR(Nlohmann json headers are required.)
fi

AC_MSG_CHECKING(nlohmann json version >= 3.7.0)
AC_TRY_CPP(
  [ #include <nlohmann/json.hpp>
    #if !(NLOHMANN_JSON_VERSION_MAJOR > 3 || \
      (NLOHMANN_JSON_VERSION_MAJOR >= 3 && NLOHMANN_JSON_VERSION_MINOR >= 7))
    #error nlohmann json too old
    #endif
  ],
  [ AC_MSG_RESULT([OK, detected nlohmann json >= 3.7.0]) ],
  [ AC_MSG_ERROR([Error: nlohman json >= 3.7.0 required])]
)

AC_SUBST(JSON_CPPFLAGS)

AC_LANG_POP(C++)

# documentation
AC_ARG_WITH(doc,
  [  --with-doc         generate documentation],
  with_doc="${withval}", with_doc="no")
AC_MSG_CHECKING(--with-doc)
DOC_SUBDIR=
DOXYGEN=
if test "$with_doc" == "yes" ; then
  AC_MSG_RESULT([yes])
  DOC_SUBDIR="doc"
  AC_CHECK_PROG(DOXYGEN, doxygen, doxygen, "")
  if ! test "$DOXYGEN" ; then
    AC_MSG_WARN(Doxygen documentation will not be built.)
  fi
else
  AC_MSG_RESULT([no])
fi
AC_SUBST(DOC_SUBDIR)
AM_CONDITIONAL([HAVE_DOXYGEN],
  [test -n "$DOXYGEN"])

# this is a cmake variable. Used in doc/doxygen/loki-b.doxy.in
PROJECT_SOURCE_DIR=`cd ${srcdir} && pwd`
AC_SUBST(PROJECT_SOURCE_DIR)

AC_OUTPUT
