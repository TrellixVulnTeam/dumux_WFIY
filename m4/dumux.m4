dnl -*- autoconf -*-
# Macros needed to find dumux and dependent libraries.  They are called by
# the macros in ${top_src_dir}/dependencies.m4, which is generated by
# "dunecontrol autogen"

# the SET_PARDISO macro moved to dune-common as of dune 2.0!
dnl AC_DEFUN([SET_PARDISO], [
dnl     AC_PREREQ(2.50)
dnl     AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
dnl     acx_pardiso_ok=no

dnl     AC_ARG_WITH(pardiso,
dnl 	[AC_HELP_STRING([--with-pardiso=<lib>], [use PARDISO library <lib>])])
dnl         case $with_pardiso in
dnl 	     yes | "") ;;
dnl              no) acx_pardiso_ok=disable ;;
dnl 	     -* | */* | *.a | *.so | *.so.* | *.o) PARDISO_LIBS="$with_pardiso" ;;
dnl 	     *) PARDISO_LIBS="-l$with_pardiso" ;;
dnl         esac

dnl     # Get fortran linker names of PARDISO functions to check for.
dnl     AC_F77_FUNC(pardisoinit)

dnl     acx_pardiso_save_LIBS="$LIBS"
dnl     LIBS="$BLAS_LIBS $LIBS $FLIBS"

dnl     # First, check PARDISO_LIBS environment variable
dnl     if test $acx_pardiso_ok = no; then
dnl     if test "x$PARDISO_LIBS" != x; then
dnl 	save_LIBS="$LIBS"; LIBS="$PARDISO_LIBS $LIBS"
dnl 	AC_MSG_CHECKING([for $pardisoinit in $PARDISO_LIBS])
dnl 	AC_TRY_LINK_FUNC($pardisoinit, [acx_pardiso_ok=yes], [PARDISO_LIBS=""])
dnl 	AC_MSG_RESULT($acx_pardiso_ok)
dnl 	LIBS="$save_LIBS"
dnl     fi
dnl     fi

dnl     AC_SUBST(PARDISO_LIBS)

dnl     LIBS="$acx_pardiso_save_LIBS"

dnl     # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
dnl     if test x"$acx_pardiso_ok" = xyes; then
dnl         ALL_PKG_LIBS="$ALL_PKG_LIBS $PARDISO_LIBS"
dnl         ifelse([$1],,AC_DEFINE(HAVE_PARDISO,1,[Define if you have a PARDISO library.]),[$1])
dnl         :
dnl     else
dnl         acx_pardiso_ok=no
dnl         $2
dnl     fi

dnl     DUNE_ADD_SUMMARY_ENTRY([Pardiso],[$acx_pardiso_ok])
dnl ]) # SET_PARDISO

# Additional checks needed to build dumux
# This macro should be invoked by every module which depends on dumux, as
# well as by dumux itself
AC_DEFUN([DUMUX_CHECKS],
[
  AC_CHECK_PROGS([TEX4HT], [tex4ht], [true])
  AC_CHECK_PROGS([MK4HT], [mk4ht], [true])
  AC_CHECK_PROGS([T4HT], [t4ht], [true])
  AC_CHECK_PROGS([LATEX], [latex], [true])
  AM_CONDITIONAL([TEX4HT], [test "x$TEX4HT" != xtrue])
  AC_CHECK_PROGS([CONVERT], [convert], [false])
  AM_CONDITIONAL([CONVERT], [test "x$CONVERT" != xfalse])

  AC_CHECK_HEADER([valgrind/memcheck.h], 
                  [HAVE_VALGRIND_H="1"],
                  AC_MSG_WARN([valgrind/memcheck.h not found]))
  AS_IF([test "x$HAVE_VALGRIND_H" = "x1"],[
    AC_DEFINE(HAVE_VALGRIND, 1, [Define whether the valgrind header files for client requests are present.])
    ])
  if test "$HAVE_VALGRIND_H" == "1"; then
     DUNE_ADD_SUMMARY_ENTRY([Valgrind client requests],["yes"])
  else
     DUNE_ADD_SUMMARY_ENTRY([Valgrind client requests],["no"])
  fi
  
  DUMUX_CHECK_QUAD

  # check whether the constexpr keyword is present
  AC_REQUIRE([CONSTEXPR_CHECK])
  # define constexpr as const if it is not available. this is quite a HACK!
  if test "x$HAVE_CONSTEXPR" != "xyes"; then
      AC_DEFINE(constexpr, const, ['set 'constexpr' to 'const' if constexpr is not supported])
  fi

  # the Boost c++ template libraries
  AX_BOOST_BASE([1.33.1])

  # check whether pardiso is installed (this is already in dune-common
  # so we only use the results.)
  # SET_PARDISO
  DUNE_ADD_SUMMARY_ENTRY([Pardiso],[$acx_pardiso_ok])

  # Add the latex and handbook status to the summary
  have_latex=no
  if test "x$LATEX" != "x"; then
     have_latex=yes
  fi
  DUNE_ADD_SUMMARY_ENTRY([Latex],[$have_latex])

  # only build the handbook if the documentation is build, latex is
  # available and the tree is checked out via a version control system
  build_handbook=no
  if test "x$enable_documentation" != "xno" && \
     test -a "$(pwd)/${top_srcdir}/stamp-vc" && \
     test "$have_latex" == "yes"; then
    build_handbook=yes
  fi
  AM_CONDITIONAL([BUILD_HANDBOOK], [test "$build_handbook" == "yes"])

  DUNE_ADD_SUMMARY_ENTRY([Build DuMuX handbook],[$build_handbook])
])

AC_DEFUN([DUMUX_CHECK_ALL],
[
DUNE_CHECK_ALL
])

AC_DEFUN([DUMUX_SUMMARY_ALL],
[
DUNE_SUMMARY_ALL
])

# DUMUX_CHECK_MODULES(NAME, HEADER, SYMBOL)
#
# THIS IS JUST A COPY OF THE DUNE_CHECK_MODULES MACRO WITH THE REQUIREMENT THAT
# THE HEADER RESIDES IN THE 'dune' SUB-DIRECTORY REMOVED!
#
# Generic check for dune modules.  This macro should not be used directly, but
# in the modules m4/{module}.m4 in the {MODULE}_CHECK_MODULE macro.  The
# {MODULE}_CHECK_MODULE macro knows the parameters to call this
# DUNE_CHECK_MODULES macro with, and it does not take any parameters itself,
# so it may be used with AC_REQUIRE.
#
# NAME   Name of the module, lowercase with dashes (like "dune-common").  The
#        value must be known when autoconf runs, so shell variables in the
#        value are not permissible.
#
# HEADER Header to check for.  The check will really be for <{HEADER}>,
#        so the header must reside within a directory called "dune".
#
# SYMBOL Symbol to check for in the module's library.  If this argument is
#        empty or missing, it is assumed that the module does not provide a
#        library.  The value must be known when autoconf runs, so shell
#        variables in the value are not permissible.  This value is actually
#        handed to AC_TRY_LINK unchanged as the FUNCTION-BODY argument, so it
#        may contain more complex stuff than a simple symbol.
#
#        The name of the library is assumed to be the same as the module name,
#        with any occurance of "-" removed.  The path of the library is
#        obtained from pkgconfig for installed modules, or assumed to be the
#        directory "lib" within the modules root for non-installed modules.
#
# In the following, {module} is {NAME} with any "-" replaced by "_" and
# {MODULE} is the uppercase version of {module}.
#
# configure options:
#   --with-{NAME}
#
# configure/shell variables:
#   {MODULE}_ROOT, {MODULE}_LIBDIR
#   HAVE_{MODULE} (1 or 0)
#   with_{module} ("yes" or "no")
#   DUNE_CPPFLAGS, DUNE_LDFLAGS, DUNE_LIBS (adds the modules values here,
#         substitution done by DUNE_CHECK_ALL)
#   ALL_PKG_CPPFLAGS, ALL_PKG_LDFLAGS, ALL_PKG_LIBS (adds the modules values
#         here, substitution done by DUNE_CHECK_ALL)
#   DUNE_PKG_CPPFLAGS, DUNE_PKG_LDFLAGS, DUNE_PKG_LIBS (deprecated, adds the
#         modules values here)
#   {MODULE}_VERSION
#   {MODULE}_VERSION_MAJOR
#   {MODULE}_VERSION_MINOR
#   {MODULE}_VERSION_REVISION
#
# configure substitutions/makefile variables:
#   {MODULE}_CPPFLAGS, {MODULE}_LDFLAGS, {MODULE}_LIBS
#   {MODULE}_ROOT
#   {MODULE}_LIBDIR (only if modules provides a library)
#
# preprocessor defines:
#   HAVE_{MODULE} (1 or undefined)
#   {MODULE}_VERSION
#   {MODULE}_VERSION_MAJOR
#   {MODULE}_VERSION_MINOR
#   {MODULE}_VERSION_REVISION
#
# automake conditionals:
#   HAVE_{MODULE}
AC_DEFUN([DUMUX_CHECK_MODULES],[
  AC_REQUIRE([AC_PROG_CXX])
  AC_REQUIRE([AC_PROG_CXXCPP])
  AC_REQUIRE([PKG_PROG_PKG_CONFIG])
  AC_REQUIRE([DUNE_DISABLE_LIBCHECK])
  AC_REQUIRE([LT_OUTPUT])

  # ____DUNE_CHECK_MODULES_____ ($1)

  m4_pushdef([_dune_name], [$1])
  m4_pushdef([_dune_module], [m4_translit(_dune_name, [-], [_])])
  m4_pushdef([_dune_header], [$2])
  m4_pushdef([_dune_ldpath], [lib])
  m4_pushdef([_dune_lib],    [m4_translit(_dune_name, [-], [])])
  m4_pushdef([_dune_symbol], [$3])
  m4_pushdef([_DUNE_MODULE], [m4_toupper(_dune_module)])

  # switch tests to c++
  AC_LANG_PUSH([C++])

  # the usual option...
  AC_ARG_WITH(_dune_name,
    AS_HELP_STRING([--with-_dune_name=PATH],[_dune_module directory]))

  # backup of flags
  ac_save_CPPFLAGS="$CPPFLAGS"
  ac_save_LIBS="$LIBS"
  ac_save_LDFLAGS="$LDFLAGS"
  CPPFLAGS=""
  LIBS=""

  ##
  ## Where is the module $1?
  ##

  AC_MSG_CHECKING([for $1 installation or source tree])

  # is a directory set?
  AS_IF([test -z "$with_[]_dune_module"],[
    #
    # search module $1 via pkg-config
    #
    with_[]_dune_module="global installation"
    AS_IF([test -z "$PKG_CONFIG"],[
      AC_MSG_RESULT([failed])
      AC_MSG_NOTICE([could not search for module _dune_name])
      AC_MSG_ERROR([pkg-config is required for using installed modules])
    ])
    AS_IF(AC_RUN_LOG([$PKG_CONFIG --exists --print-errors "$1"]),[
      _dune_cm_CPPFLAGS="`$PKG_CONFIG --cflags _dune_name`" 2>/dev/null
      _DUNE_MODULE[]_ROOT="`$PKG_CONFIG --variable=prefix _dune_name`" 2>/dev/null 
      _DUNE_MODULE[]_VERSION="`$PKG_CONFIG --modversion _dune_name`" 2>/dev/null
      _dune_cm_LDFLAGS=""
      ifelse(_dune_symbol,,
        [_DUNE_MODULE[]_LIBDIR=""
         _dune_cm_LIBS=""],
        [_DUNE_MODULE[]_LIBDIR=`$PKG_CONFIG --variable=libdir _dune_name 2>/dev/null`
         _dune_cm_LIBS="-L$_DUNE_MODULE[]_LIBDIR -l[]_dune_lib"])
      HAVE_[]_DUNE_MODULE=1
      AC_MSG_RESULT([global installation in $_DUNE_MODULE[]_ROOT])
    ],[
      HAVE_[]_DUNE_MODULE=0
      AC_MSG_RESULT([not found])
    ])
  ],[
    #
    # path for module $1 is specified via command line
    #
    AS_IF([test -d "$with_[]_dune_module"],[
      # expand tilde / other stuff
      _DUNE_MODULE[]_ROOT=`cd $with_[]_dune_module && pwd`

      # expand search path (otherwise empty CPPFLAGS)
      AS_IF([test -d "$_DUNE_MODULE[]_ROOT/include/dune"],[
        # Dune was installed into directory given by with-dunecommon
        _dune_cm_CPPFLAGS="-I$_DUNE_MODULE[]_ROOT/include"
        _DUNE_MODULE[]_BUILDDIR=_DUNE_MODULE[]_ROOT
        _DUNE_MODULE[]_VERSION="`PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$_DUNE_MODULE[]_ROOT/lib/pkgconfig $PKG_CONFIG --modversion _dune_name`" 2>/dev/null
      ],[
        _DUNE_MODULE[]_SRCDIR=$_DUNE_MODULE[]_ROOT
        # extract src and build path from Makefile, if found
	    AS_IF([test -f $_DUNE_MODULE[]_ROOT/Makefile],[
          _DUNE_MODULE[]_SRCDIR="`sed -ne '/^abs_top_srcdir = /{s/^abs_top_srcdir = //; p;}' $_DUNE_MODULE[]_ROOT/Makefile`"
		])
        _dune_cm_CPPFLAGS="-I$_DUNE_MODULE[]_SRCDIR"
        _DUNE_MODULE[]_VERSION="`grep Version $_DUNE_MODULE[]_SRCDIR/dune.module | sed -e 's/^Version: *//'`" 2>/dev/null
      ])
      _dune_cm_LDFLAGS=""
      ifelse(_dune_symbol,,
        [_DUNE_MODULE[]_LIBDIR=""
         _dune_cm_LIBS=""],
        [_DUNE_MODULE[]_LIBDIR="$_DUNE_MODULE[]_ROOT/lib"
         _dune_cm_LIBS="-L$_DUNE_MODULE[]_LIBDIR -l[]_dune_lib"])
      # set expanded module path
      with_[]_dune_module="$_DUNE_MODULE[]_ROOT"
      HAVE_[]_DUNE_MODULE=1
      AC_MSG_RESULT([found in $_DUNE_MODULE[]_ROOT])
    ],[
      HAVE_[]_DUNE_MODULE=0
      AC_MSG_RESULT([not found])
      AC_MSG_ERROR([_dune_name-directory $with_[]_dune_module does not exist])
    ])
  ])

  CPPFLAGS="$ac_save_CPPFLAGS $DUNE_CPPFLAGS $_dune_cm_CPPFLAGS"
  ##  
  ## check for an arbitrary header
  ##
  AC_CHECK_HEADER([_dune_header],
    [HAVE_[]_DUNE_MODULE=1],
    [HAVE_[]_DUNE_MODULE=0
     AS_IF([test -n "$_DUNE_MODULE[]_ROOT"],[
       AC_MSG_WARN([$_DUNE_MODULE[]_ROOT does not seem to contain a valid _dune_name (_dune_header not found)])
     ])
    ]
  )

  ##
  ## check for lib (if lib name was provided)
  ##
  ifelse(_dune_symbol,,
    AC_MSG_NOTICE([_dune_name does not provide libs]),

    AS_IF([test "x$enable_dunelibcheck" = "xno"],[
      AC_MSG_WARN([library check for _dune_name is disabled. DANGEROUS!])
    ],[
      AS_IF([test "x$HAVE_[]_DUNE_MODULE" = "x1"],[

        # save current LDFLAGS
        ac_save_CXX="$CXX"
        HAVE_[]_DUNE_MODULE=0

        # define LTCXXLINK like it will be defined in the Makefile
        CXX="./libtool --tag=CXX --mode=link $ac_save_CXX"

        # use module LDFLAGS
        LDFLAGS="$ac_save_LDFLAGS $DUNE_LDFLAGS $DUNE_PKG_LDFLAGS $_dune_cm_LDFLAGS"
        LIBS="$_dune_cm_LIBS $DUNE_LIBS $LIBS"

        AC_MSG_CHECKING([for lib[]_dune_lib])

        AC_TRY_LINK(dnl
          [#]include<_dune_header>,
          _dune_symbol,
          [
            AC_MSG_RESULT([yes])
            HAVE_[]_DUNE_MODULE=1
          ],[
            AC_MSG_RESULT([no])
            HAVE_[]_DUNE_MODULE=0
            AS_IF([test -n "$_DUNE_MODULE[]_ROOT"],[
             AC_MSG_WARN([$with_[]_dune_module does not seem to contain a valid _dune_name (failed to link with lib[]_dune_lib[].la)])
            ])
          ]
        )
      ])

      # reset variables
      CXX="$ac_save_CXX"
    ])
  )

  # did we succeed?
  AS_IF([test "x$HAVE_[]_DUNE_MODULE" = "x1"],[
    # add the module's own flags and libs to the modules and the global
    # variables
    DUNE_ADD_MODULE_DEPS(m4_defn([_dune_name]), m4_defn([_dune_name]),
        [$_dune_cm_CPPFLAGS], [$_dune_cm_LDFLAGS], [$_dune_cm_LIBS])

    # set variables for our modules
    AC_SUBST(_DUNE_MODULE[]_CPPFLAGS, "$_DUNE_MODULE[]_CPPFLAGS")
    AC_SUBST(_DUNE_MODULE[]_LDFLAGS, "$_DUNE_MODULE[]_LDFLAGS")
    AC_SUBST(_DUNE_MODULE[]_LIBS, "$_DUNE_MODULE[]_LIBS")
    AC_SUBST(_DUNE_MODULE[]_ROOT, "$_DUNE_MODULE[]_ROOT")
    ifelse(m4_defn([_dune_symbol]),,
      [],
      [AC_SUBST(_DUNE_MODULE[]_LIBDIR)
    ])
    AC_DEFINE(HAVE_[]_DUNE_MODULE, 1, [Define to 1 if] _dune_name [was found])

    DUNE_PARSE_MODULE_VERSION(_dune_name, $_DUNE_MODULE[]_VERSION)

    # set DUNE_* variables
    # This should actually be unneccesary, but I'm keeping it in here for now
    # for backward compatibility
    DUNE_LDFLAGS="$DUNE_LDFLAGS $_DUNE_MODULE[]_LDFLAGS"
    DUNE_LIBS="$_DUNE_MODULE[]_LIBS $DUNE_LIBS"
    
    # add to global list
    # only add my flags other flags are added by other packages 
    DUNE_PKG_CPPFLAGS="$DUNE_PKG_CPPFLAGS $_DUNE_MODULE[]_CPPFLAGS"
    DUNE_PKG_LIBS="$DUNE_PKG_LIBS $LIBS"
    DUNE_PKG_LDFLAGS="$DUNE_PKG_LDFLAGS $_DUNE_MODULE[]_LDFLAGS"

    with_[]_dune_module="yes"
  ],[
    with_[]_dune_module="no"
  ])

  AM_CONDITIONAL(HAVE_[]_DUNE_MODULE, test x$HAVE_[]_DUNE_MODULE = x1)

  # reset previous flags
  CPPFLAGS="$ac_save_CPPFLAGS"
  LDFLAGS="$ac_save_LDFLAGS"
  LIBS="$ac_save_LIBS"

  # add this module to DUNE_SUMMARY
  DUNE_MODULE_ADD_SUMMARY_ENTRY(_dune_name)

  # remove local variables
  m4_popdef([_dune_name])
  m4_popdef([_dune_module])
  m4_popdef([_dune_header])
  m4_popdef([_dune_ldpath])
  m4_popdef([_dune_lib])
  m4_popdef([_dune_symbol])
  m4_popdef([_DUNE_MODULE])

  # restore previous language settings (leave C++)
  AC_LANG_POP([C++])
])

# Additional checks needed to find dumux
# This macro should be invoked by every module which depends on dumux, but
# not by dumux itself
AC_DEFUN([DUMUX_CHECK_MODULE],
[
  DUMUX_CHECK_MODULES([dumux],[dumux/common/exceptions.hh])
])
