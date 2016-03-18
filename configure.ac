#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#
# NEED TO DISABLE THE OPTION:  --disable-dependency-tracking
# did not work without it on Farideh's mac. Bug of automake. 
#
AC_PREREQ([2.68])
LT_PREREQ([2.4])
AC_INIT([biogeme], [2.4], [http://groups.yahoo.com/group/biogeme],[biogeme],[http://biogeme.epfl.ch])
AC_LANG(C++)
AC_CONFIG_AUX_DIR([.])
AC_CONFIG_MACRO_DIR([m4])
LT_INIT([disable-static])
#LT_INIT
LT_LANG([C++])

AM_INIT_AUTOMAKE([no-dependencies -Wno-portability])
AC_CONFIG_FILES([Makefile 
	libraries/utils/Makefile 
	libraries/linearAlgebra/Makefile 
	libraries/trustRegion/Makefile 
	libraries/cfsqp/Makefile 
	libraries/donlp2/Makefile 
	libraries/solvopt/Makefile 
	libraries/randomNumbers/Makefile 
	libraries/gaussHermite/Makefile 
	libraries/parameters/Makefile 
	libraries/bisonbiogeme/Makefile 
	libraries/pythonbiogeme/Makefile 
	libraries/Makefile 
        tools/histograms/Makefile
        tools/prepareData/Makefile
        tools/mod2py/Makefile
        biogeme/Makefile
	biosim/Makefile 
	biomerge/Makefile 
        pythonbiogeme/Makefile        
	pybiogeme/Makefile
	testing/Makefile
        doc/Makefile
        webui/Makefile
        gui/Makefile])
AC_CONFIG_FILES([pythonbiogeme/pythonbiogeme],[chmod +x pythonbiogeme/pythonbiogeme])
AC_CONFIG_FILES([pythonbiogeme/debugbiogeme],[chmod +x pythonbiogeme/debugbiogeme])
AC_CONFIG_FILES([webui/webbiogeme],[chmod +x webui/webbiogeme])
AC_CONFIG_FILES([tools/prepareData/biocheckdata],[chmod +x tools/prepareData/biocheckdata])
AC_CONFIG_FILES([tools/prepareData/biosampledata],[chmod +x tools/prepareData/biosampledata])
AC_CONFIG_FILES([tools/prepareData/biopreparedata],[chmod +x tools/prepareData/biopreparedata])
AC_CONFIG_FILES([tools/histograms/histogram],[chmod +x tools/histograms/histogram])
AC_CONFIG_FILES([tools/histograms/weightedhistogram],[chmod +x tools/histograms/weightedhistogram])
AC_CONFIG_FILES([tools/mod2py/mod2py],[chmod +x tools/mod2py/mod2py])
AC_CONFIG_FILES([pythonbiogeme/pythonbiogeme_static],[chmod +x pythonbiogeme/pythonbiogeme_static])
AC_CONFIG_SRCDIR([libraries/utils/patAbsTime.cc])
AC_CONFIG_HEADERS([config.h])

##########################################################################
# debug compilation support
##########################################################################

# Support for extended precision
AC_ARG_ENABLE([longdouble],AC_HELP_STRING([--enable-longdouble],[Use long double floating point data type (slower)])
)

AC_MSG_CHECKING([whether to build with debug information])
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
        [enable debug data generation (def=no)])],
    [debugit="$enableval"],
    [debugit=no])
AC_MSG_RESULT([$debugit])

if test x"$debugit" = x"yes"; then
    AC_DEFINE([DEBUG],[],[Debug Mode])
#    AM_CXXFLAGS="$AM_CXXFLAGS -g -Wall -Werror -Wno-uninitialized -Wno-unused-variable -O0"
    AM_CXXFLAGS="$AM_CXXFLAGS -g -Wall -O0 -fno-omit-frame-pointer"
    AM_LDFLAGS="$AM_LDFLAGS -g"
    WITHDEBUG="yes"
else
    AC_DEFINE([NODEBUG],[],[No-debug Mode])
    AM_CXXFLAGS="$AM_CXXFLAGS -O3"
    WITHDEBUG="no (use option --enable-debug)"
fi
AC_SUBST([AM_CXXFLAGS])

##########################################################################


#Check OS. This is needed to see if yywrap must be added or not. There
#must be a better way, but I did not find it.

# Check for running on Darwin.
#UNAME=`uname -s`
#AC_MSG_CHECKING([if running on Darwin])
#if test "${UNAME}"="Darwin"; then
#  AC_DEFINE([IS_DARWIN], 1, [running on Darwin])
#  AC_MSG_RESULT(yes)
#  AC_SUBST([libsuffix],[".dylib"])
#  AM_CONDITIONAL([DARWIN], true)
#else
#  AC_MSG_RESULT(no)
#  AC_SUBST([libsuffix],[".so"])
#  AM_CONDITIONAL([DARWIN], false)
#fi
#AC_MSG_CHECKING([if running on Mingw])
#if test "${UNAME}"="MINGW32_NT-5.1"; then
#  AC_DEFINE([IS_MINGW], 1, [running on Mingw])
#  AC_MSG_RESULT(yes)
#  AM_CONDITIONAL([MINGW], true)
#else
#  AC_MSG_RESULT(no)
#  AM_CONDITIONAL([MINGW], false)
#fi

AC_ARG_ENABLE([gui],
  AS_HELP_STRING([--enable-gui],[Build the version with a GUI. Requires fltk.])
)
AC_ARG_ENABLE([python],
  AS_HELP_STRING([--enable-python],[Build the python version. Requires python3.1, python3.2, python3.3, python3.3m, python3.4 or python3.5.])
)
AC_ARG_ENABLE([bison],
  AS_HELP_STRING([--enable-bison],[Build the regular version.])
)

AC_ARG_ENABLE([staticexec],
  AS_HELP_STRING([--enable-staticexec],[Build executables with static libraries.])
)

# Support for versioned installs
AC_ARG_ENABLE([versioned],AC_HELP_STRING([--enable-versioned],[Install with versioned filenames and paths])
)


#AC_ARG_VAR(HOSTTYPE,(Apple only) Select one: i386 x86_64 ppc ppc64)

# Check if GNU scientific library is available

#AC_CHECK_LIB([m],[cos])
#AC_CHECK_LIB([gslcblas],[cblas_dgemm])
#AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

#if test "x$ac_cv_lib_gsl_gsl_blas_dgemm" == "xyes"; then 
#     WITHGSL="yes"
#else
#     WITHGSL="no"
#fi


# custom program locations
# python31
AC_ARG_WITH([python_dir],
  AS_HELP_STRING([--with-python-dir=DIR],[Specify the directory where Python 3.1, 3.2, 3.3, 3.3m, 3.4 or 3.5 is.]),
  [python_path=$withval],
  [python_path=$PATH]
)

# fltk
AC_ARG_WITH([fltk_dir],
  AS_HELP_STRING([--with-fltk-dir=DIR],[Specify the directory where fltk is installed.]),
  [fltk_path=$withval],
  [fltk_path=$PATH]
)

# bison++
AC_ARG_WITH([bison_dir],
  AS_HELP_STRING([--with-bison-dir=DIR],[Specify the directory where bison++ is installed.]),
  [bison_path=$withval],
  [bison_path=$PATH]
)

# lldb
AC_ARG_WITH([lldb_dir],
  AS_HELP_STRING([--with-lldb-dir=DIR],[Specify the directory where lldb is installed.]),
  [lldb_path=$withval],
  [lldb_path=$PATH]
)

# flex++
AC_ARG_WITH([flex_dir],
  AS_HELP_STRING([--with-flex-dir=DIR],[Specify the directory where flex++ is installed.]),
  [flex_path=$withval],
  [flex_path=$PATH]
)

# latexmk
AC_ARG_WITH([latexmk_dir],
  AS_HELP_STRING([--with-latexmk-dir=DIR],[Specify the directory where latexmk is installed.]),
  [latexmk_path=$withval],
  [latexmk_path=$PATH]
)


AC_PROG_CXX([clang++ g++])
AC_PROG_CC([clang gcc])
# Remove the next statement before distribution
AC_DISABLE_SHARED
AC_PROG_LIBTOOL

AC_PATH_PROG([FLTKCONFIG],[fltk-config],[no],[$fltk_path])
if test "x$FLTKCONFIG" = "xno"; then
  WINBIO="no"
  if test "x$with_fltk_dir" = "x" ; then 
    WINBIO_SUGGEST="(Download and install fltk)"
  else
    WINBIO_SUGGEST="(fltk-config not found in $with_fltk_dir)"
  fi
#  FLTKCONFIG="$SHELL $missing_dir/missing fltk-config"
else
  if test "x$enable_gui" = "xyes" ; then
    WINBIO="yes"
    WINBIO_SUGGEST=""
    FLTK_LDFLAGS="`$FLTKCONFIG --ldflags`"
    FLTK_LDSTATICFLAGS="`$FLTKCONFIG --ldstaticflags`"
    FLTK_CXXFLAGS="`$FLTKCONFIG --cxxflags`"
    AC_SUBST(FLTK_LDFLAGS)
    AC_SUBST(FLTK_LDSTATICFLAGS)
    AC_SUBST(FLTK_CXXFLAGS)
  else
    WINBIO="no"
    WINBIO_SUGGEST="(use option --enable-gui)"
  fi
fi

AC_PATH_PROGS([PYTHONBIN],[python3.5 python3.4 python3.3m python3.3 python3.2 python3.1],[no],[$python_path])
if test "x$PYTHONBIN" = "xno"; then
  NEWBIO="no"
  if test "x$with_python_dir" = "x" ; then
    NEWBIO_SUGGEST="(Download and install python3.5, python3.4, python3.3m, python 3.3, python 3.2 or python 3.1)"
  else
    NEWBIO_SUGGEST="(python 3.1, 3.2, 3.3, 3.3, 3.4 or 3.5 not found in $with_python_dir)"
  fi
#  PYTHONBIN = "$SHELL $missing_dir/missing python3.2"
else 
  PYTHONBINBASE="`basename $PYTHONBIN`"
  case "${PYTHONBINBASE}" in
   "python3.5") PYTHONCONFIG="$srcdir/pythonbiogeme/python/biogeme-python3.5-config.py";;
   "python3.4") PYTHONCONFIG="$srcdir/pythonbiogeme/python/biogeme-python3.4-config.py";;
   "python3.3m") PYTHONCONFIG="$srcdir/pythonbiogeme/python/biogeme-python3.3-config.py";;
   "python3.3") PYTHONCONFIG="$srcdir/pythonbiogeme/python/biogeme-python3.3-config.py";;
   "python3.2") PYTHONCONFIG="$srcdir/pythonbiogeme/python/biogeme-python3.2-config.py";;
   "python3.1") PYTHONCONFIG="$srcdir/pythonbiogeme/python/biogeme-python3.1-config.py";;
   "python3.2-cygwin") PYTHONCONFIG="$srcdir/pythonbiogeme/python/biogeme-python3.2-cygwin-config.py";;
  esac
  PYTHONHEADERDIR="`$PYTHONBIN $PYTHONCONFIG headerdir`"
  AC_CHECK_HEADERS([$PYTHONHEADERDIR/Python.h],[python_header=yes],[python_header=no])
  if test "x$enable_python" = "xyes" -a "x$python_header" = "xyes"; then
    NEWBIO="yes"
    NEWBIO_SUGGEST="(Python version `$PYTHONBIN $PYTHONCONFIG version`)"
    PYTHON_INCLUDES="`$PYTHONBIN $PYTHONCONFIG includes`"
    PYTHON_CXXFLAGS="`$PYTHONBIN $PYTHONCONFIG cxxflags`"
    PYTHON_LIBS="`$PYTHONBIN $PYTHONCONFIG libs`"
    PYTHON_LDFLAGS="`$PYTHONBIN $PYTHONCONFIG ldflags`"
    AC_SUBST(PYTHON_INCLUDES)
    AC_SUBST(PYTHON_CXXFLAGS)
    AC_SUBST(PYTHON_LIBS)
    AC_SUBST(PYTHON_LDFLAGS)
    AC_SUBST(PYTHONBIN)
  else 
    NEWBIO="no"
    if test "x$python_header" = "xno"; then
      NEWBIO_SUGGEST="(Install python3-dev)"
    else
      NEWBIO_SUGGEST="(use option --enable-python)"
    fi
  fi
fi

if test "x$enable_bison" = "xyes" ; then
  BISONBIO="yes"
  BISONBIO_SUGGEST=""
else 
  BISONBIO="no"
  BISONBIO_SUGGEST="(use option --enable-bison)"
fi

AC_PATH_PROG([BISON],[bison++],[bison++],[$bison_path])
AM_MISSING_PROG([YACC], [$BISON])

AC_PATH_PROG([LLDB],[lldb],[lldb],[$lldb_path])

AC_PATH_PROG([FLEX],[flex++],[flex++],[$flex_path])
AM_MISSING_PROG([LEX], [$FLEX])

AC_PATH_PROG([LATEXMK],[latexmk],[:],[$latexmk_path])
AC_SUBST([LEX_OUTPUT_ROOT],[lex.yy])



# Support for portable batch system

AC_PATH_PROG([PBS],[qsub],[no])
AM_CONDITIONAL([SYSTEM_WITH_PBS],[test "x$PBS" != "xno"])

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h sys/time.h errno.h math.h malloc.h sys/malloc.h stdio.h pthread.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_FUNC_MKTIME
AC_CHECK_FUNCS([floor ceil gettimeofday localtime_r localtime time pow exp log sqrt atan mktime])
AC_CHECK_FUNC(chdir,
  [AC_DEFINE([HAVE_CHDIR], [1], [Define to 1 if you have chdir.])],
  [AC_CHECK_FUNC(SetCurrentDirectory,
    [AC_DEFINE([HAVE_SETCURRENTDIRECTORY], [1], [Define to 1 if you have SetCurrentDirectory.])])]
)

#Decide what to build

AM_CONDITIONAL([BUILD_DEBUG],[test "x$debugit" = "xyes"])
AM_CONDITIONAL([BUILD_STATIC], [test "x$enable_staticexec" = "xyes"])
AM_CONDITIONAL([BUILD_LONGDOUBLE], [test "x$enable_longdouble" = "xyes"])
AM_CONDITIONAL([BUILD_GUI],[test "x$FLTKCONFIG" != "xno" -a "x$enable_gui" = "xyes"])
AM_CONDITIONAL([BUILD_PYTHON],[test "x$PYTHONBIN" != "xno" -a "x$enable_python" = "xyes"])
AM_CONDITIONAL([BUILD_BISONBIOGEME],[test "x$enable_bison" = "xyes"])

if test "x$enable_staticexec" = "xyes" ; then
  STATIC="yes"
else
  STATIC="no (use option --enable-staticexec)"
fi

if test "x$enable_longdouble" = "xyes" ; then
  WITHLONGDOUBLE="yes"
  AC_DEFINE([LONGDOUBLE],[1],[Use long double floating point data type])
else
  WITHLONGDOUBLE="no (use option --enable-longdouble)"
fi

# Check if it is an Apple computer. This is to insert Apple specific options
AM_CONDITIONAL([APPLE], [AS_CASE([$host_os], [darwin*], [true], [false])])

#Decide where to install

AM_CONDITIONAL([INSTALL_VERSIONED], [test "x$enable_versioned" = "xyes" ])
if test "x$enable_versioned" = "xyes" ; then
  INSTALL_VERSION_SUFFIX="-${VERSION}"
  INSTALL_RELEASE="-release ${VERSION}"
  VERSIONED="yes (with suffix -${VERSION})"
else
  INSTALL_VERSION_SUFFIX=""
  INSTALL_RELEASE=""
  VERSIONED="no (use option --enable-versioned)"
fi
AC_SUBST(INSTALL_RELEASE)
AC_SUBST([INSTALL_VERSION_SUFFIX])

AC_OUTPUT

# GNU Scientific library...$WITHGSL

AC_MSG_RESULT([
------------------------------------------------------------------------
  $PACKAGE $VERSION:  Automatic configuration OK.
 Build biogeme............$BISONBIO $BISONBIO_SUGGEST
 Build biomerge...........$BISONBIO $BISONBIO_SUGGEST
 Build biosim.............$BISONBIO $BISONBIO_SUGGEST
 Build python biogeme.....$NEWBIO $NEWBIO_SUGGEST
 Build gui version........$WINBIO $WINBIO_SUGGEST
 Static executables.......$STATIC
 Long double..............$WITHLONGDOUBLE
 Executables for debug....$WITHDEBUG
 Versioned installation...$VERSIONED
])