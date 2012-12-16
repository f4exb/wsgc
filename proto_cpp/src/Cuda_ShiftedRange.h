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

 shifted iterator (or range) - inspired by the strided range

 */

#ifndef __CUDA_SHIFTED_RANGE_H__
#define __CUDA_SHIFTED_RANGE_H__

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/functional.h>

/**
 * \brief shifted values in a range
 *
 * These examples illustrate how to make shifted values access within a range of values:
 *   shifted_range([0, 1, 2, 3, 4, 5, 6], 0)  -> [0, 1, 2, 3, 4, 5, 6]
 *   shifted_range([0, 1, 2, 3, 4, 5, 6], 1)  -> [1, 2, 3, 4, 5, 6, 0]
 *   shifted_range([0, 1, 2, 3, 4, 5, 6], 2)  -> [2, 3, 4, 5, 6, 0, 1]
 *   shifted_range([0, 1, 2, 3, 4, 5, 6], -1) -> [6, 0, 1, 2, 3, 4, 5]
 *   ...
 */
template<typename Iterator>
class shifted_range
{
public:

	typedef typename thrust::iterator_difference<Iterator>::type difference_type;

	struct shift_functor: public thrust::unary_function<difference_type,
			difference_type>
	{
		difference_type size;
		difference_type shift;

		shift_functor(difference_type size, difference_type shift) :
			size(size), shift(shift)
		{
		}

		__host__ __device__
		difference_type operator()(const difference_type& i) const
		{
			return (i + shift + size) % size;
		}
	};

	typedef typename thrust::counting_iterator<difference_type> CountingIterator;
	typedef typename thrust::transform_iterator<shift_functor, CountingIterator> TransformIterator;
	typedef typename thrust::permutation_iterator<Iterator, TransformIterator> PermutationIterator;

	// type of the repeat_range iterator
	typedef PermutationIterator iterator;

	// construct repeat_range for the range [first,last)
	shifted_range(Iterator first, Iterator last, difference_type shift) :
			first(first), last(last), shift(shift)
	{
	}

	iterator begin(void) const
	{
		return PermutationIterator(first,
				TransformIterator(CountingIterator(0), shift_functor((last - first), shift)));
	}

	iterator end(void) const
	{
		return begin() + (last - first);
	}

protected:
	Iterator first;
	Iterator last;
	difference_type shift;
};

#endif // __CUDA_SHIFTED_RANGE_H__
