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

    Collection of index transforms

*/

#ifndef __CUDA_INDEX_TRANSFORMS_H__
#define __CUDA_INDEX_TRANSFORMS_H__

#include <cutil_inline.h>    // includes cuda.h and cuda_runtime_api.h


/**
 * \brief Convert a linear index ranging [0..T*Isf*Ihf[ to index for source FFT
*/
struct transpose_index_A : public thrust::unary_function<size_t,size_t>
{
	size_t _T;   //!< FFT size
	size_t _Isf; //!< Number of frequency sub-steps (a.k.a. sub-frequencies)

	__host__ __device__
	transpose_index_A() : _T(4096), _Isf(8) {}

	__host__ __device__
	transpose_index_A(size_t T, size_t Isf) : _T(T), _Isf(Isf) {}

	__host__ __device__
	size_t operator()(size_t linear_index)
	{
		size_t t = linear_index % _T;
		size_t isf = (linear_index / _T) % _Isf;
		return (_T * isf) + t;
	}
};


/**
 * \brief Convert a linear index ranging [0..T*Isf*Ihf[ to index for PRN local code FFT conjugate
*/
struct transpose_index_B : public thrust::unary_function<size_t,size_t>
{
	size_t _T;   //!< FFT size
	size_t _Isf; //!< Number of frequency sub-steps or sub-frequencies
	size_t _Ihf; //!< Number of frequency steps or harmonic frequencies

	__host__ __device__
	transpose_index_B() : _T(4096), _Isf(8), _Ihf(127) {}

	__host__ __device__
	transpose_index_B(size_t T, size_t Isf, size_t Ihf) : _T(T), _Isf(Isf), _Ihf(Ihf) {}

	__host__ __device__
	size_t operator()(size_t linear_index)
	{
		unsigned int fhi = linear_index / (_Isf * _T);
		unsigned int t = linear_index % _T;
		return (t + fhi - (_Ihf/2)) % _T;
	}
};


/**
 * \brief Convert a linear index ranging [0..T*Isf*Ihf[ to index for IFFT input storage
*/
struct transpose_index_C : public thrust::unary_function<size_t,size_t>
{
	size_t _B;   //!< Batch size
	size_t _b;   //!< Batch index

	__host__ __device__
	transpose_index_C() : _B(3), _b(0) {}

	__host__ __device__
	transpose_index_C(size_t B, size_t b) : _B(B), _b(b) {}

	__host__ __device__
	size_t operator()(size_t linear_index)
	{
		return _b + linear_index*2*_B;
	}
};


/**
 * \brief Convert a linear index ranging [0..T*Isf*Ihf*B[ to index for IFFT output averaging
 */
struct transpose_index_avgC : public thrust::unary_function<size_t,size_t>
{
	size_t _B;       //!< Batch size (half buffer size)
	size_t _T;       //!< FFT size
	size_t _Isf;     //!< Number of frequency sub-steps or sub-frequencies
	size_t _Ihf;     //!< Number of frequency steps or harmonic frequencies
	size_t _ais;     //!< Average index at start of batch sequence

	__host__ __device__
	transpose_index_avgC() : _B(3), _T(4096), _Isf(8), _Ihf(127), _ais(0) {}

	__host__ __device__
	transpose_index_avgC(size_t B, size_t T, size_t Isf, size_t Ihf, size_t ais, size_t b_shift) : _B(B), _T(T), _Isf(Isf), _Ihf(Ihf), _ais(ais) {}

	__host__ __device__
	size_t operator()(size_t linear_index)
	{
		/*
		size_t t = linear_index % _T;
		size_t isf = (linear_index / _T) % _Isf;
		size_t ihf = (linear_index / _T) / _Isf;
		unsigned int z = t + _T*isf + _T*_Isf*ihf; // zero average index position
		unsigned int ai = (_ais + _b_shift + (linear_index % _B)) % (2*_B); // current average index
		return z + ai;
		*/
		//return ((linear_index % _T) + _T*((linear_index / _T) % _Isf) + _T*_Isf*((linear_index / _T) / _Isf)) + ((_ais + _b_shift + (linear_index % _B)) % (2*_B));
		return ((_ais + (linear_index / (_T*_Isf*_Ihf))) % (2*_B)) + 2*_B*(linear_index % (_T*_Isf*_Ihf));
	}
};


/**
 * \brief Convert a linear index ranging [0..T*Isf*Ihf*B[ to index for IFFT output averaging - DOES NOT WORK!
 * What it is supposed to do:
 * _B=3 _ais=0;
 * |0|1|2|3|..|6| 7| 8| 9|......|1|2|3|4|..| 7| 8| 9|10| < one series corresponds to one complete average calculation
 * |1|2|3|4|..|7| 8| 9|10|......|2|3|4|5|..| 8| 9|10|11| | But it works this way
 * |2|3|4|5|..|8| 9|10|11|......|3|4|5|0|..| 9|10|11| 6| | => the last one is reused!
 * |3|4|5|0|..|9|10|11| 6|......|4|5|0|1|..|10|11| 6| 7| V
 */
struct transpose_index_multi_stride : public thrust::unary_function<size_t,size_t>
{
	size_t _B;       //!< Batch size (half buffer size)
	size_t _T;       //!< FFT size
	size_t _Isf;     //!< Number of frequency sub-steps or sub-frequencies
	size_t _Ihf;     //!< Number of frequency steps or harmonic frequencies
	size_t _ais;     //!< Average index at start of batch sequence

	__host__ __device__
	transpose_index_multi_stride() : _B(3), _T(4096), _Isf(8), _Ihf(127), _ais(0) {}

	__host__ __device__
	transpose_index_multi_stride(size_t B, size_t T, size_t Isf, size_t Ihf, size_t ais) : _B(B), _T(T), _Isf(Isf), _Ihf(Ihf), _ais(ais) {}

	__host__ __device__
	size_t operator()(size_t linear_index)
	{
		return ((_ais + (linear_index / (_T*_Isf*_Ihf))) % (2*_B)) + 2*_B*(linear_index % (_T*_Isf*_Ihf));
	}
};


/**
 * \brief Convert a linear index ranging [0..T*Isf*Ihf*B[ to index for IFFT output averaging. So called "skipping" stride.
 * This one works for average calculation
 * Example with _B = 3, _ais = 0:
 * |0|1|2|...|6|7|8|...|12|13|14|... < one series corresponds to one member of the average calculation (now you only have _B groups sizes)
 *        ^
 *        skips
 * with _ais = 5 (last "odd" items in the batch):
 * |5|0|1|...|11|6|7|...|15|12|13|...
 *
 * _B is the batch size so the stride of FFT is actually 2*_B. The "odd" batches are skipped.
 * "odd" batches are done when moving the start on an odd batch, then the even are skipped.
 * _ais is the origin shift. It corresponds to the index in a 2 batches set
 */
struct transpose_index_skipping_stride : public thrust::unary_function<size_t,size_t>
{
	size_t _B;       //!< Batch size (half buffer size)
	size_t _ais;     //!< Average index at start of batch sequence

	__host__ __device__
	transpose_index_skipping_stride() : _B(3), _ais(0) {}

	__host__ __device__
	transpose_index_skipping_stride(size_t B, size_t ais) : _B(B), _ais(ais) {}

	__host__ __device__
	size_t operator()(size_t linear_index)
	{
		//return ((_ais + linear_index) % (2*_B)) + ((linear_index / (2*_B))*2*_B);
        return ((_ais + (linear_index % _B)) % (2*_B)) + ((linear_index/_B)*2*_B);
	}
};


#endif // __CUDA_INDEX_TRANSFORMS_H__
