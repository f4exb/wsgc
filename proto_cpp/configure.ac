AC_INIT(wsgc_test, 0.1.0)
AC_CONFIG_MACRO_DIR([m4])
m4_include([m4/boost.m4])
AM_INIT_AUTOMAKE(wsgc_test, 0.1.0)
AC_LANG_PUSH([C++])

ARCH=`uname -m`

# ------------------------------------------------------------------------------
# Boost
# ------------------------------------------------------------------------------
BOOST_REQUIRE([1.42])

# ------------------------------------------------------------------------------
# FFTW3
# ------------------------------------------------------------------------------

FFTW3_REQUIRED_VERSION=3.3.2
PKG_CHECK_MODULES(FFTW3F, fftw3f >= $FFTW3_REQUIRED_VERSION)
AC_SUBST(FFTW3F_CFLAGS)
AC_SUBST(FFTW3F_LIBS)

# ------------------------------------------------------------------------------
# DOXYGEN
# ------------------------------------------------------------------------------

AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN";
   then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi

AM_CONDITIONAL([HAVE_DOXYGEN], 
[test -n "$DOXYGEN"])AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([doc/Doxyfile])])

# ------------------------------------------------------------------------------
# RSSOFT
# ------------------------------------------------------------------------------

# Setup defaut RSSoft paths
RSSOFT_INSTALL_PATH="/usr/local/rssoft"

# Check: enable RSSoft?
AC_MSG_CHECKING([whether to use RSSoft library])
AC_ARG_ENABLE(rssoft,
    [AS_HELP_STRING(--enable-rssoft,enable RSSoft support)],
    rssoftval="yes",
    rssoftval="no")
AC_MSG_RESULT([$rssoftval])

# check to see if RSSoft been selected
if test "x$rssoftval" == "xno"; then
    AC_MSG_WARN([Not using RSSoft!])
fi

# Setup custom RSSoft paths
AC_ARG_WITH([rssoft],
   [AS_HELP_STRING(
        [--with-rssoft=PATH],
        [prefix where RSSoft is installed @<:@default=auto@:>@])],
   [RSSOFT_INSTALL_PATH=$withval],
   [with_rssoft=auto])

# Verify RSSoft paths
if test  "x$rssoftval" == "xyes"; then
    AC_CHECK_FILE([${RSSOFT_INSTALL_PATH}/include/GFq.h],[hasrssoft="yes"],[hasrssoft="no"])
    if test "x$hasrssoft" == "xno"; then
       AC_MSG_WARN([RSSoft library not found - will not have RSSoft support])
       rssoftval="no"
    fi
    AM_CONDITIONAL([HAS_RSSOFT],[test "x$rssoftval" == "xyes"]) 
else
    AM_CONDITIONAL([HAS_RSSOFT],[false])     
fi

# Setup RSSoft flags
if test  "x$rssoftval" == "xyes"; then
	RSSOFT_CFLAGS="-I$RSSOFT_INSTALL_PATH/include"
	RSSOFT_LDFLAGS="-L$RSSOFT_INSTALL_PATH/lib -lrssoft"
	AC_SUBST(RSSOFT_CFLAGS)
	AC_SUBST(RSSOFT_LDFLAGS)
	AC_DEFINE([_RSSOFT],[1],[Defined if RSSoft should be used])
fi


# ------------------------------------------------------------------------------
# CCSOFT
# ------------------------------------------------------------------------------

# Setup defaut CCSoft paths
RSSOFT_INSTALL_PATH="/usr/local/ccsoft"

# Check: enable CCSoft?
AC_MSG_CHECKING([whether to use CCSoft library])
AC_ARG_ENABLE(ccsoft,
    [AS_HELP_STRING(--enable-ccsoft,enable CCSoft support)],
    ccsoftval="yes",
    ccsoftval="no")
AC_MSG_RESULT([$ccsoftval])

# check to see if CCSoft been selected
if test "x$ccsoftval" == "xno"; then
    AC_MSG_WARN([Not using CCSoft!])
fi

# Setup custom RSSoft paths
AC_ARG_WITH([ccsoft],
   [AS_HELP_STRING(
        [--with-ccsoft=PATH],
        [prefix where CCSoft is installed @<:@default=auto@:>@])],
   [CCSOFT_INSTALL_PATH=$withval],
   [with_ccsoft=auto])

# Verify CCSoft paths
if test  "x$ccsoftval" == "xyes"; then
    AC_CHECK_FILE([${CCSOFT_INSTALL_PATH}/include/CC_SequentialDecoding.h],[hasccsoft="yes"],[hasccsoft="no"])
    if test "x$hasccsoft" == "xno"; then
       AC_MSG_WARN([CCSoft library not found - will not have CCSoft support])
       ccsoftval="no"
    fi
    AM_CONDITIONAL([HAS_CCSOFT],[test "x$ccsoftval" == "xyes"]) 
else
    AM_CONDITIONAL([HAS_CCSOFT],[false])     
fi

# Setup RSSoft flags
if test  "x$ccsoftval" == "xyes"; then
    CCSOFT_CFLAGS="-I$CCSOFT_INSTALL_PATH/include"
    CCSOFT_LDFLAGS="-L$CCSOFT_INSTALL_PATH/lib -lccsoft"
    AC_SUBST(CCSOFT_CFLAGS)
    AC_SUBST(CCSOFT_LDFLAGS)
    AC_DEFINE([_CCSOFT],[1],[Defined if CCSoft should be used])
fi


# ------------------------------------------------------------------------------
# CUDA
# ------------------------------------------------------------------------------

# Setup default CUDA paths
CUDA_INSTALL_PATH="/usr/local/cuda"
CUDA_SDK_INSTALL_PATH="/usr/local/NVIDIA_GPU_Computing_SDK"
CUDA_THRUST_INSTALL_PATH="/usr/local/cuda/thrust"

# Check: enable CUDA?
AC_MSG_CHECKING([whether to use CUDA library])
AC_ARG_ENABLE(cuda,
    [AS_HELP_STRING(--enable-cuda,enable CUDA support)],
    cudaval="yes",
    cudaval="no")
AC_MSG_RESULT([$cudaval])

# check to see if CUDA been selected
if test "x$cudaval" == "xno"; then
    AC_MSG_WARN([Not using CUDA!])
fi

# Setup custom CUDA paths
AC_ARG_WITH([cuda],
   [AS_HELP_STRING(
        [--with-cuda=PATH],
        [prefix where CUDA is installed @<:@default=auto@:>@])],
   [CUDA_INSTALL_PATH=$withval],
   [with_cuda=auto])
AC_ARG_WITH([cuda-sdk],
   [AS_HELP_STRING(
        [--with-cuda-sdk=PATH],
        [prefix where CUDA SDK is installed @<:@default=auto@:>@])],
   [CUDA_SDK_INSTALL_PATH=$withval],
   [with_cuda_sdk=auto])
   
# Setup optional target hardware architecture
AC_ARG_WITH([cuda-arch],
   [AS_HELP_STRING(
        [--with-cuda-arch=ARCH],
        [CUDA hardware architecture passed as -arch= argument to NVCC compiler @<:@default=none@:>@])],
   [CUDA_ARCH="-arch=$withval"],
   [with_cuda_arch=""])

# Check: enable CUDA Emulator?
AC_ARG_ENABLE([cudaemu],[--enable-cuda-emu Turn on device emulation for CUDA],
    [case "${enableval}" in
        yes) cudaemulation=true;;
        no)  cudaemulation=false;;
        *) AC_MSG_ERROR([bad value ${enableval} for --enable-cuda-emu]);;
    esac],
    [cudaemulation=false])


# Verify CUDA paths

if test  "x$cudaval" == "xyes"; then
    AC_CHECK_PROGS([hasnvcc], [nvcc], [nonvcc], "$CUDA_INSTALL_PATH/bin")
    if test "x$hasnvcc" == "xnonvcc"; then
       AC_MSG_WARN([nvcc not found - will not have CUDA support])
       cudaval="no"
    fi;
    AC_CHECK_FILE([${CUDA_SDK_INSTALL_PATH}/C/common/inc/cutil.h],[hascudasdk="yes"],[hascudasdk="no"])
    if test "x$hascudasdk" == "xno"; then
       AC_MSG_WARN([CUDA SDK not found - will not have CUDA support])
       cudaval="no"
    fi    
    AC_CHECK_FILE([${CUDA_INSTALL_PATH}/include/thrust/tuple.h],[hascudathrust="yes"],[hascudathrust="no"])
    if test "x$hascudathrust" == "xno"; then
       AC_MSG_WARN([CUDA Thrust not found - will not have CUDA support])
       cudaval="no"
    fi
    AM_CONDITIONAL([HAS_CUDA],[test "x$cudaval" == "xyes"])     
else
    AM_CONDITIONAL([HAS_CUDA],[false])     
fi
    
# Setup nvcc flags

NVCCFLAGS="$NVCCFLAGS -ccbin gcc-4.5"

case "$ARCH" in
    amd64*)
        CUDA_LIB_SUFFIX="lib64"
    ;;
    x86_64*)
        CUDA_LIB_SUFFIX="lib64"
    ;;
    *)
        CUDA_LIB_SUFFIX=""
    ;;
esac

if test "x$cudaemulation" == "xtrue"; then
    NVCCFLAGS="$NVCCFLAGS -deviceemu"
fi

# initialize array index
N="0"

if test  "x$cudaval" == "xyes"; then
    NVCCFLAGS="$NVCCFLAGS"
    CUDA_CFLAGS="$CUDA_CFLAGS"
    CUDA_CFLAGS="$CUDA_CFLAGS -I$CUDA_SDK_INSTALL_PATH/C/common/inc/"
    CUDA_CFLAGS="$CUDA_CFLAGS -I$CUDA_SDK_INSTALL_PATH/shared/inc/"
    CUDA_CFLAGS="$CUDA_CFLAGS -I$CUDA_INSTALL_PATH/include"
    CUDA_LDFLAGS="-L$CUDA_INSTALL_PATH/$CUDA_LIB_SUFFIX"
    CUDA_LIBS="-lcuda -lcudart -lcufft -lcublas"
    NVCCFLAGS="$NVCCFLAGS $CUDA_CFLAGS"
    NVLDFLAGS="$CUDA_LDFLAGS $CUDA_LIBS"
    NVCC="$CUDA_INSTALL_PATH/bin/nvcc"
	if test  "x$rssoftval" == "xyes"; then
		NVCCFLAGS="$NVCCFLAGS -D_RSSOFT"
	fi
    AC_SUBST(CUDA_CFLAGS)
    AC_SUBST(NVCC)
    AC_SUBST(NVCCFLAGS)
    AC_SUBST(NVLDFLAGS)
    AC_DEFINE([_CUDA],[1],[Defined if CUDA should be used])
    eval "ARRAY${N}='-D_CUDA'"
    N=`expr $N + 1`
fi


# construct compiler flags
N=`expr $N - 1`
TMP0=
eval TMP0="\$ARRAY0"
for i in $(seq 1 $N)
do
    eval TMP0="\$TMP0' '\$ARRAY${i}"
done

# set C compiler flags
AC_SUBST(CDFLAGS, $TMP0)
    
# ------------------------------------------------------------------------------
# EPILOG
# ------------------------------------------------------------------------------
    
dnl: this is a comment
dnl: ${CXXFLAGS=-g}
AC_PROG_CXX

AC_CONFIG_FILES([Makefile doc/Makefile src/Makefile])
AC_OUTPUT


