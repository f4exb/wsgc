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

 Repeat with increment iterator (or range) - inspired by thrust example strided_range.cu

 */

#ifndef __CUDA_REPEAT_INCREMENTAL_RANGE_H__
#define __CUDA_REPEAT_INCREMENTAL_RANGE_H__

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/functional.h>

/**
 * \brief Repeated incremental range or iterator
 *
 * This example illustrates how to make repeated incremental access to a range of values
 *
 * Examples:
 *   repeat_incremental_range([0, 1, 2, 3, 4, 5, 6], 4) -> [0, 1, 2, 3,  1, 2, 3, 4,  2, 3, 4, 5,  3, 4, 5, 6]
 *   ...
 */
template<typename Iterator>
class repeat_incremental_range
{
public:

	typedef typename thrust::iterator_difference<Iterator>::type difference_type;

	struct repeat_incremental_functor: public thrust::unary_function<difference_type,
			difference_type>
	{
		difference_type repeat;

		repeat_incremental_functor(difference_type repeat) :
				repeat(repeat)
		{
		}

		__host__ __device__
		difference_type operator()(const difference_type& i) const
		{
            return (i%repeat) + (i/repeat);
		}
	};

	typedef typename thrust::counting_iterator<difference_type> CountingIterator;
	typedef typename thrust::transform_iterator<repeat_incremental_functor, CountingIterator> TransformIterator;
	typedef typename thrust::permutation_iterator<Iterator, TransformIterator> PermutationIterator;

	// type of the repeat_incremental_range iterator
	typedef PermutationIterator iterator;

	// construct repeat_incremental_range for the range [first,last)
	repeat_incremental_range(Iterator first, Iterator last, difference_type repeat) :
			first(first), last(last), repeat(repeat)
	{
	}

	iterator begin(void) const
	{
		return PermutationIterator(first,
				TransformIterator(CountingIterator(0), repeat_incremental_functor(repeat)));
	}

	iterator end(void) const
	{
		return begin() + ((last - first) - repeat + 1) * repeat; 
	}

protected:
	Iterator first;
	Iterator last;
	difference_type repeat;
};

#endif // __CUDA_REPEAT_INCREMENTAL_RANGE_H__
