dnl libgd checker for autoconf.
dnl
dnl Cheks for the gd/png/zlib library and all necessary stuff.
dnl

AC_DEFUN([GD_CHECK],
#
# Handle user hints
#
[
 AC_MSG_CHECKING([if libgd is wanted])
 
 AC_ARG_WITH([libgd],
  [  --with-libgd=DIR root directory path of libgd installation
		   defaults to /usr and /usr/local
  --without-libgd to disable libgd usage completely (not yet implemented)],
  [
    #
    # Run this if -with or -without was specified
    #
    if test "$withval" != no ; then
       AC_MSG_RESULT(yes)
       LIBGD_WANTED=yes
       if test "$withval" != yes ; then
         LIBGD_HOME_DIR="$withval"
       fi
    else
       LIBGD_WANTED=no
       AC_MSG_RESULT(no)
    fi]
  ,[
   # Nothing was said I assume libgd is needed
   
    LIBGD_WANTED=yes
    AC_MSG_RESULT(yes)
  ])



# just do the job if libgd is wanted
if test $LIBGD_WANTED = yes; then
 
    # save the state at begining
    _old_LDFLAGS=$LDFLAGS
    _old_CPPFLAGS=$CPPFLAGS

    # perform search in various directory.
    AC_MSG_CHECKING([for libgd with png support])
    for _search_dir_to_inc in "$LIBGD_HOME_DIR"  /usr /usr/local ; do

	LDFLAGS="$_old_LDFLAGS -L${_search_dir_to_inc}/lib"
	CPPFLAGS="$_old_CPPFLAGS -I${_search_dir_to_inc}/include"

    
	AC_TRY_COMPILE([#include "gd.h"] , 
		       [ gdImagePtr im;
			 FILE *pngout;
			 gdImagePng(im, pngout);
		       ] ,
		       _compil_ok="yes" )  
	if test "$_compil_ok" = yes ; then
	    break
	fi
    done

    # if OK, then just add the library to $LIBS, 
    # else reset to the initial state

    if test "$_compil_ok" = yes ; then
	AC_MSG_RESULT([yes])
	AC_DEFINE(HAVE_LIBGD, 1, is libgd present)
	LIBS="$LIBS -lgd -lpng -lz"
	if  test -z "$_search_dir_to_inc" ; then
       	   LDFLAGS=$_old_LDFLAGS
	   CPPFLAGS=$_old_CPPFLAGS
	fi
    else
	AC_MSG_RESULT([no])
	AC_MSG_WARN([libgd not found, graphic topologies will be unavailable])
    fi
    AM_CONDITIONAL(USE_GD_LIB_SRC, test "$_compil_ok")
    # if libgd is not wanted, disable some piece of code
    # still to implement
    # I'm thinking about 
    # setting up a #define and check for it in the code.

fi
])




