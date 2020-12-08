AC_PREREQ([2.69])
AC_INIT([loki-b],
	[m4_esyscmd_s([cat VERSION])],
	[loki@tecnico.ulisboa.pt])
AC_CONFIG_SRCDIR([Makefile.am])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([autodir])
AC_CONFIG_MACRO_DIR([autodir])
AM_INIT_AUTOMAKE([1.11 foreign subdir-objects])
AC_CONFIG_FILES([
	Makefile
	lib/Makefile
	doc/Makefile
	doc/doxygen/Makefile
	doc/doxygen/Doxyfile
	testsuite/Makefile
	])

LT_INIT([disable-static win32-dll])
AC_SUBST(LIBTOOL_DEPS)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for header files.

# eigen
AC_LANG_PUSH(C++)
AC_ARG_WITH(eigen,
    [  --with-eigen(=eigen prefix)
              use the Eigen library],
    with_eigen="${withval}", with_eigen="yes")
AC_MSG_CHECKING(--with-eigen)
AC_MSG_RESULT($with_eigen)
EIGEN_CPPFLAGS=
if test "$with_eigen" != "no" ; then
  dnl Testsuites that rely on eigen
  if test "$with_eigen" = "yes" ; then
    AC_CHECK_HEADER(eigen3/Eigen/Core,
      with_eigen=yes,
      AC_MSG_ERROR(The eigen3 headers were not found. Please install these and reconfigure.))
  else
    EIGEN_CPPFLAGS="-I$with_eigen"
    CPPFLAGS="$CPPFLAGS $EIGEN_CPPFLAGS"
    AC_CHECK_HEADER(eigen3/Eigen/Core,
    with_eigen=yes,
    AC_MSG_ERROR(The eigen3 headers were not found. Please install these and reconfigure.))
  fi
  AC_DEFINE_UNQUOTED([HAVE_EIGEN],1,[Define if we have the Eigen matrix vector library])
else
  AC_MSG_ERROR(Eigen3 headers are required.)
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

# set CMAKE variable
CMAKE_CURRENT_BINARY_DIR=`pwd`
AC_SUBST(CMAKE_CURRENT_BINARY_DIR)
CMAKE_CURRENT_SOURCE_DIR=`cd ${srcdir} && pwd`
AC_SUBST(CMAKE_CURRENT_SOURCE_DIR)

AC_OUTPUT
