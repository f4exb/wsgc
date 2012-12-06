AC_INIT(wsgc_test, 0.1.0)
AM_INIT_AUTOMAKE(wsgc_test, 0.1.0)
AC_LANG_PUSH([C++])

ARCH=`uname -m`

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


