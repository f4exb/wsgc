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

 Repeating iterator (or range) - inspired by the strided range
 This is different of repeated_range example (see repeat_value for this)
 Here the whole range repeats itself like the tiled_range but with only one
 differnce_type argument

 */

#ifndef __CUDA_REPEAT_RANGE_H__
#define __CUDA_REPEAT_RANGE_H__

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/functional.h>

/**
 * \brief Repetition of the range of values
 *
 * These examples illustrate how to make repeated access to a range of values
 *   repeat_range([0, 1, 2], 1) -> [0, 1, 2]
 *   repeat_range([0, 1, 2], 2) -> [0, 1, 2, 0, 1, 2]
 *   repeat_range([0, 1, 2], 3) -> [0, 1, 2, 0, 1, 2, 0, 1, 2]
 *   ...
 */
template<typename Iterator>
class repeat_range
{
public:

	typedef typename thrust::iterator_difference<Iterator>::type difference_type;

	struct modulo_functor: public thrust::unary_function<difference_type,
			difference_type>
	{
		difference_type modulo;

		modulo_functor(difference_type modulo) :
			modulo(modulo)
		{
		}

		__host__ __device__
		difference_type operator()(const difference_type& i) const
		{
			return i % modulo;
		}
	};

	typedef typename thrust::counting_iterator<difference_type> CountingIterator;
	typedef typename thrust::transform_iterator<modulo_functor, CountingIterator> TransformIterator;
	typedef typename thrust::permutation_iterator<Iterator, TransformIterator> PermutationIterator;

	// type of the repeat_range iterator
	typedef PermutationIterator iterator;

	// construct repeat_range for the range [first,last)
	repeat_range(Iterator first, Iterator last, difference_type repeat) :
			first(first), last(last), repeat(repeat)
	{
	}

	iterator begin(void) const
	{
		return PermutationIterator(first,
				TransformIterator(CountingIterator(0), modulo_functor(last - first)));
	}

	iterator end(void) const
	{
		return begin() + ((last - first) * repeat);
	}

protected:
	Iterator first;
	Iterator last;
	difference_type repeat;
};

#endif // __CUDA_REPEAT_RANGE_H__
