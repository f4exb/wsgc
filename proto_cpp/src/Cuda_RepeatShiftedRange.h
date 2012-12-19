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

 Repeat with shift iterator (or range) - inspired by thrust example strided_range.cu

 */

#ifndef __CUDA_REPEAT_SHIFTED_RANGE_H__
#define __CUDA_REPEAT_SHIFTED_RANGE_H__

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/functional.h>

/**
 * \brief Repeated shifted range or iterator
 *
 * This example illustrates how to make repeated shifted access to a range of values
 *
 * Examples:
 *   repeat_shifted_range([0, 1, 2, 3, 4, 5, 6], -1, 3) -> [6, 0, 1, 2, 3, 4, 5,  0, 1, 2, 3, 4, 5, 6,  1, 2, 3, 4, 5, 6, 0 ] : initial shift -1, repetition 3
 *   ...
 */
template<typename Iterator>
class repeat_shifted_range
{
public:

	typedef typename thrust::iterator_difference<Iterator>::type difference_type;

	struct repeat_shifted_functor: public thrust::unary_function<difference_type,
			difference_type>
	{
		difference_type length;
		difference_type init_shift;

		repeat_shifted_functor(difference_type length, difference_type init_shift) :
				length(length), init_shift(init_shift)
		{
		}

		__host__ __device__
		difference_type operator()(const difference_type& i) const
		{
            return (((i%length) + (i/length)) + init_shift + length) % length;
		}
	};

	typedef typename thrust::counting_iterator<difference_type> CountingIterator;
	typedef typename thrust::transform_iterator<repeat_shifted_functor, CountingIterator> TransformIterator;
	typedef typename thrust::permutation_iterator<Iterator, TransformIterator> PermutationIterator;

	// type of the repeat_shifted_range iterator
	typedef PermutationIterator iterator;

	// construct repeat_shifted_range for the range [first,last)
	repeat_shifted_range(Iterator first, Iterator last, difference_type init_shift, difference_type repeat) :
			first(first), last(last), init_shift(init_shift), repeat(repeat)
	{
	}

	iterator begin(void) const
	{
		return PermutationIterator(first,
				TransformIterator(CountingIterator(0), repeat_shifted_functor((last-first), init_shift)));
	}

	iterator end(void) const
	{
		return begin() + ((last - first) * repeat); 
	}

protected:
	Iterator first;
	Iterator last;
	difference_type init_shift;
	difference_type repeat;
};

#endif // __CUDA_REPEAT_SHIFTED_RANGE_H__
