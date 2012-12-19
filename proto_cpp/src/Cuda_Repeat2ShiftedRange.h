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

 Repeating repeated (repeat-repeat or repeat squared) with shift iterator (or range) - inspired by thrust example strided_range.cu

 */

#ifndef __CUDA_REPEAT_2_SHIFTED_RANGE_H__
#define __CUDA_REPEAT_2_SHIFTED_RANGE_H__

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/functional.h>

/**
 * \brief Repeating repeated shifted range or iterator
 *
 * This example illustrates how to make repeating repeated shifted access to a range of values
 *
 * Example initial shift -1, sub-repetition (repeating same) 2, repetition 3
 *   repeat_shifted_range([0, 1, 2, 3, 4, 5, 6], -1, 2, 3) -> 
 *        [6, 0, 1, 2, 3, 4, 5,  6, 0, 1, 2, 3, 4, 5,  
 *         0, 1, 2, 3, 4, 5, 6,  0, 1, 2, 3, 4, 5, 6,
 *         1, 2, 3, 4, 5, 6, 0,  1, 2, 3, 4, 5, 6, 0 ] 
 *   ...
 */
template<typename Iterator>
class repeat_2_shifted_range
{
public:

	typedef typename thrust::iterator_difference<Iterator>::type difference_type;

	struct repeat_2_shifted_functor: public thrust::unary_function<difference_type,
			difference_type>
	{
		difference_type length;
		difference_type init_shift;
        difference_type sub_repeat;

		repeat_2_shifted_functor(difference_type length, difference_type init_shift, difference_type sub_repeat) :
				length(length), init_shift(init_shift), sub_repeat(sub_repeat)
		{
		}

		__host__ __device__
		difference_type operator()(const difference_type& i) const
		{
            return ((i%length) + (i/(length*sub_repeat)) + init_shift + length) % length;
		}
	};

	typedef typename thrust::counting_iterator<difference_type> CountingIterator;
	typedef typename thrust::transform_iterator<repeat_2_shifted_functor, CountingIterator> TransformIterator;
	typedef typename thrust::permutation_iterator<Iterator, TransformIterator> PermutationIterator;

	// type of the repeat_2_shifted_range iterator
	typedef PermutationIterator iterator;

	// construct repeat_2_shifted_range for the range [first,last)
	repeat_2_shifted_range(Iterator first, Iterator last, difference_type init_shift, difference_type sub_repeat, difference_type repeat) :
			first(first), last(last), init_shift(init_shift), sub_repeat(sub_repeat), repeat(repeat)
	{
	}

	iterator begin(void) const
	{
		return PermutationIterator(first,
				TransformIterator(CountingIterator(0), repeat_2_shifted_functor((last-first), init_shift, sub_repeat)));
	}

	iterator end(void) const
	{
		return begin() + ((last - first) * sub_repeat * repeat); 
	}

protected:
	Iterator first;
	Iterator last;
	difference_type init_shift;
	difference_type sub_repeat;
	difference_type repeat;
};

#endif // __CUDA_REPEAT_2_SHIFTED_RANGE_H__
