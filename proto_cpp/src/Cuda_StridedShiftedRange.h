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

 Strided iterator (or range) - found in Thrust library examples
 http://code.google.com/p/thrust/source/browse/examples/strided_range.cu

 */

#ifndef __CUDA_STRIDED_SHIFTED_RANGE_H__
#define __CUDA_STRIDED_SHIFTED_RANGE_H__

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/functional.h>

/**
 * \brief Strided and shifted range of values
 *
 * These examples illustrate how to make strided and shifted access to a range of values:
 *   strided_shifted_range([0, 1, 2, 3, 4, 5, 6, 7 ,8], 3, 0) -> [0, 3, 6] : this is like the strided_range with stride 3
 *   strided_shifted_range([0, 1, 2, 3, 4, 5, 6, 7 ,8], 3, 1) -> [1, 4, 7]
 *   strided_shifted_range([0, 1, 2, 3, 4, 5, 6, 7 ,8], 3, 2) -> [2, 5, 8]
 */
template<typename Iterator>
class strided_shifted_range
{
public:

	typedef typename thrust::iterator_difference<Iterator>::type difference_type;

	struct stride_shift_functor: public thrust::unary_function<difference_type,
			difference_type>
	{
		difference_type stride;
		difference_type shift;

		stride_shift_functor(difference_type stride, difference_type shift) :
				stride(stride), shift(shift)
		{
		}

		__host__ __device__
		difference_type operator()(const difference_type& i) const
		{
			return stride * i + shift;
		}
	};

	typedef typename thrust::counting_iterator<difference_type> CountingIterator;
	typedef typename thrust::transform_iterator<stride_shift_functor, CountingIterator> TransformIterator;
	typedef typename thrust::permutation_iterator<Iterator, TransformIterator> PermutationIterator;

	// type of the strided_shifted_range iterator
	typedef PermutationIterator iterator;

	// construct strided_shifted_range for the range [first,last)
	strided_shifted_range(Iterator first, Iterator last, difference_type stride, difference_type shift) :
			first(first), last(last), stride(stride), shift(shift)
	{
	}

	iterator begin(void) const
	{
		return PermutationIterator(first,
				TransformIterator(CountingIterator(0), stride_shift_functor(stride ,shift)));
	}

	iterator end(void) const
	{
		return begin() + ((last - first) + (stride - 1)) / stride;
	}

protected:
	Iterator first;
	Iterator last;
	difference_type stride;
	difference_type shift;
};

#endif // __CUDA_STRIDED_SHIFTED_RANGE_H__
