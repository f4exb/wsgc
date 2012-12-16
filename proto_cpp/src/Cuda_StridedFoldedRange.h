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

 Strided folded iterator (or range) - inspired by thrust example strided_range.cu

 */

#ifndef __CUDA_STRIDED_FOLDED_RANGE_H__
#define __CUDA_STRIDED_FOLDED_RANGE_H__

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/functional.h>

/**
 * \brief Strided folded range or iterator
 *
 * This example illustrates how to make strided folded access to a range of values
 *
 * Examples:
 *   strided_folded_range([0, 1, 2, 3, 4, 5, 6], 1) -> [0, 1, 2, 3, 4, 5, 6]
 *   strided_folded_range([0, 1, 2, 3, 4, 5, 6], 2) -> [0, 2, 4, 6, 1, 3, 5]
 *   strided_folded_range([0, 1, 2, 3, 4, 5, 6], 3) -> [0, 3, 6, 1, 4, 2, 5]
 *   strided_folded_range([0, 1, 2, 3, 4, 5, 6, 7, 8], 3) -> [0, 3, 6, 1, 4, 7, 2, 5, 8]
 *   strided_folded_range([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], 4) -> [0, 4, 8, 1, 5, 9, 2, 6, 10, 3, 7, 11]
 *   ...
 * As you can see in these examples the size of the range is reduced to the nearest multiple of the stride that is
 * smaller than the original size. This is because we don't want any logic in the process just a calculation of the
 * new index for best performance with CUDA Thrust. The result of the last example on the compete range would be:
 * [0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11] and this cannot be done without checking if the strided value
 * would be larger than the size.
 */
template<typename Iterator>
class strided_folded_range
{
public:

	typedef typename thrust::iterator_difference<Iterator>::type difference_type;

	struct stride_fold_functor: public thrust::unary_function<difference_type,
			difference_type>
	{
		difference_type stride;
		difference_type usable_size;

		stride_fold_functor(difference_type stride, difference_type usable_size) :
				stride(stride), usable_size(usable_size)
		{
		}

		__host__ __device__
		difference_type operator()(const difference_type& i) const
		{
            return (i*stride)%usable_size + (i*stride)/usable_size ;
		}
	};

	typedef typename thrust::counting_iterator<difference_type> CountingIterator;
	typedef typename thrust::transform_iterator<stride_fold_functor, CountingIterator> TransformIterator;
	typedef typename thrust::permutation_iterator<Iterator, TransformIterator> PermutationIterator;

	// type of the strided_range iterator
	typedef PermutationIterator iterator;

	// construct strided_range for the range [first,last)
	strided_folded_range(Iterator first, Iterator last, difference_type stride) :
			first(first), last(last), stride(stride), usable_size(((last - first)/stride)*stride)
	{
	}

	iterator begin(void) const
	{
		return PermutationIterator(first,
				TransformIterator(CountingIterator(0), stride_fold_functor(stride, usable_size)));
	}

	iterator end(void) const
	{
		return begin() + usable_size; 
	}

protected:
	Iterator first;
	Iterator last;
	difference_type stride;
	difference_type usable_size;
};

#endif // __CUDA_STRIDED_FOLDED_RANGE_H__
