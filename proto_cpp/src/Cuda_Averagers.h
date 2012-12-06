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

    Averaging specialized classes

*/

#ifndef __CUDA_AVERAGERS_H__
#define __CUDA_AVERAGERS_H__

#include "Cuda_IndexTransforms.h"
#include "Cuda_Operators.h"

typedef struct AveragingDimensions_s
{
	size_t _B;    //!< Batch size
	size_t _T;    //!< FFT/IFFT size
	size_t _Ifs;  //!< Number of sub-frequency steps or sub-frequencies
	size_t _Ifh;  //!< Number of frequency steps or harmonic frequencies
} AveragingDimensions_t;

template <size_t N>
class Averager
{
public:
	void run(thrust::device_vector<cuComplex>& mP, AveragingDimensions_t& a, size_t b_shift)
	{}
};

/*
 * \brief Averager on 4 values
 */
/*
template <>
class Averager<4>
{
public:
	void run(thrust::device_vector<cuComplex>& mP, AveragingDimensions_t& a, size_t b_shift)
	{
		thrust::for_each(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 1+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 2+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 3+b_shift)))
				)
			),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 1+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 2+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 3+b_shift)))
				)
			),
			avgsum_functor_complex<4>()
			//null_operator()
		);
	}
};
*/

template <>
class Averager<4>
{
public:
	void run(thrust::device_vector<cuComplex>& mP, AveragingDimensions_t& a, size_t b_shift)
	{
		thrust::for_each(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_skipping_stride(a._B, b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_skipping_stride(a._B, 1+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_skipping_stride(a._B, 2+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_skipping_stride(a._B, 3+b_shift)))
				)
			),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_skipping_stride(a._B, b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_skipping_stride(a._B, 1+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_skipping_stride(a._B, 2+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_skipping_stride(a._B, 3+b_shift)))
				)
			),
			avgsum_functor_complex<4>()
			//null_operator()
		);
	}
};


/*
 * \brief Averager on 6 values
 */
template <>
class Averager<6>
{
public:
	void run(thrust::device_vector<cuComplex>& mP, AveragingDimensions_t& a, size_t b_shift)
	{
		thrust::for_each(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 1+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 2+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 3+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 4+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 5+b_shift)))
				)
			),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 1+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 2+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 3+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 4+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 5+b_shift)))
				)
			),
			avgsum_functor_complex<6>()
		);
	}
};


/*
 * \brief Averager on 8 values
 */
template <>
class Averager<8>
{
public:
	void run(thrust::device_vector<cuComplex>& mP, AveragingDimensions_t& a, size_t b_shift)
	{
		thrust::for_each(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 1+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 2+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 3+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 4+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 5+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 6+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 7+b_shift)))
				)
			),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 1+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 2+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 3+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 4+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 5+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 6+b_shift))),
					thrust::make_permutation_iterator(mP.begin(), thrust::make_transform_iterator(thrust::make_counting_iterator(0)+(a._T*a._Ifs*a._Ifh*a._B), transpose_index_multi_stride(a._B, a._T, a._Ifs, a._Ifh, 7+b_shift)))
				)
			),
			avgsum_functor_complex<8>()
		);
	}
};


#endif // __CUDA_AVERAGERS_H__


