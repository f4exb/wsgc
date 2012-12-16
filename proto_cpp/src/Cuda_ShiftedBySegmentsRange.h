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

 shifted by (regular) segments iterator (or range) - inspired by the strided range

 */

#ifndef __CUDA_SHIFTED_BY_SEGMENTS_RANGE_H__
#define __CUDA_SHIFTED_BY_SEGMENTS_RANGE_H__

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/functional.h>

/**
 * \brief Shifted range by equal length segments:
 *
 * These examples illustrate how to make shifted access by segments of equal length:
 *   strided_shifted_range([0, 1, 2, 10, 11, 12, 20, 21, 22], 3, 1)   -> [1, 2, 0, 11, 12, 10, 21, 22, 20]
 *   strided_shifted_range([0, 1, 2, 10, 11, 12, 20, 21, 22], 3, -1)  -> [2, 0, 1, 12, 10, 11, 22, 20, 21]
 *   ...
 */
template<typename Iterator>
class shifted_by_segments_range
{
public:

	typedef typename thrust::iterator_difference<Iterator>::type difference_type;

	struct shift_by_segments_functor: public thrust::unary_function<difference_type,
			difference_type>
	{
		difference_type segment_size;
		difference_type shift;

		shift_by_segments_functor(difference_type segment_size, difference_type shift) :
			segment_size(segment_size), shift(shift)
		{
		}

		__host__ __device__
		difference_type operator()(const difference_type& i) const
		{
			return ((i/segment_size)*segment_size) + (((i%segment_size)+shift+segment_size)%segment_size);
		}
	};

	typedef typename thrust::counting_iterator<difference_type> CountingIterator;
	typedef typename thrust::transform_iterator<shift_by_segments_functor, CountingIterator> TransformIterator;
	typedef typename thrust::permutation_iterator<Iterator, TransformIterator> PermutationIterator;

	// type of the shifted_by_segments_range iterator
	typedef PermutationIterator iterator;

	// construct shifted_by_segments_range for the range [first,last)
	shifted_by_segments_range(Iterator first, Iterator last, difference_type segment_size, difference_type shift) :
			first(first), last(last), segment_size(segment_size), shift(shift)
	{
	}

	iterator begin(void) const
	{
		return PermutationIterator(first,
				TransformIterator(CountingIterator(0), shift_by_segments_functor(segment_size, shift)));
	}

	iterator end(void) const
	{
		return begin() + (last - first);
	}

protected:
	Iterator first;
	Iterator last;
	difference_type segment_size;
	difference_type shift;
};

#endif // __CUDA_SHIFTED_BY_SEGMENTS_RANGE_H__
