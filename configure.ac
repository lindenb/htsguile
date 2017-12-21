dnl If anyone wants to improve on this go right ahead
AC_INIT([htsguile],0.1)
AC_PROG_CC
AC_PROG_CC_C99
AC_CONFIG_HEADERS([config.h])

dnl Where is htslib
AC_ARG_WITH([htslib_include],
            [AC_HELP_STRING([--with-htslib-include=DIR],
                            [HTSlib include directory])],
            [CPPFLAGS="-I$withval $CPPFLAGS"])

AC_ARG_WITH([htslib_lib],
            [AC_HELP_STRING([--with-htslib-lib=DIR],
                            [HTSlib library directory])],
            [LDFLAGS="-L$withval $CPPFLAGS"])

AC_ARG_WITH([htslib],
            [AC_HELP_STRING([--with-htslib=DIR],
                            [HTSlib installation directory])],
            [
            if test -z "$with_htslib_lib" -a -z "$with_htslib_include" ; then
                CPPFLAGS="-I$withval/include $CPPFLAGS"
                LDFLAGS="-L$withval/lib $LDFLAGS"
            else
                AC_MSG_FAILURE([Do not use --with-htslib and --with-htslib-include/--with-htslib-lib options simultaneously.])
            fi
            ])

dnl Check htslib
AC_CHECK_HEADER([htslib/hts.h],
    [AC_SEARCH_LIBS([hts_open], [hts],
        [],
        [AC_MSG_ERROR([Failed to find HTSLib library])])],
    [AC_MSG_ERROR([Failed to find HTSLib headers])])


AC_OUTPUT([Makefile])
