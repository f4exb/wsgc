/*
     Copyright 2012 Edouard Griffiths <f4exb at free dot fr>
 
     This file is part of WSGC. A Weak Signal transmission mode using Gold Codes
 
     This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
 
     You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
     Foundation, Inc., 51 Franklin Street, Boston, MA  02110-1301  USA

     Static not real time prototype in C++
      
     WsgcTypes

     Collection of typedefs and defines mainly used to switch between single and double precision
*/
#ifndef __WSGC_TYPES_H__
#define __WSGC_TYPES_H__

#include <fftw3.h>
#include <complex>

#ifdef WSGC_FLOAT_DOUBLE
typedef std::complex<double> wsgc_complex;
typedef fftw_complex wsgc_fftw_complex;
typedef fftw_plan wsgc_fftw_plan;
typedef double  wsgc_float;
#define WSGC_FFTW_PLAN fftw_plan_dft_1d
#define WSGC_FFTW_PLAN_MANY fftw_plan_many_dft
#define WSGC_FFTW_DESTROY_PLAN fftw_destroy_plan
#define WSGC_FFTW_MALLOC fftw_malloc
#define WSGC_FFTW_FREE fftw_free
#define WSGC_FFTW_EXECUTE fftw_execute

#else
typedef std::complex<float> wsgc_complex;
typedef fftwf_complex wsgc_fftw_complex;
typedef fftwf_plan wsgc_fftw_plan;
typedef float  wsgc_float;
#define WSGC_FFTW_PLAN fftwf_plan_dft_1d
#define WSGC_FFTW_PLAN_MANY fftwf_plan_many_dft
#define WSGC_FFTW_DESTROY_PLAN fftwf_destroy_plan
#define WSGC_FFTW_MALLOC fftwf_malloc
#define WSGC_FFTW_FREE fftwf_free
#define WSGC_FFTW_EXECUTE fftwf_execute
#endif

#endif // __WSGC_TYPES_H__
