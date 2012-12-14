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

 Repeating values iterator (or range) - inspired by the strided range

 */

#ifndef __CUDA_REPEAT_VALUES_H__
#define __CUDA_REPEAT_VALUES_H__

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/functional.h>

#include <thrust/fill.h>
#include <thrust/device_vector.h>

// for printing
#include <thrust/copy.h>
#include <ostream>

// this example illustrates how to make repeated values with a range of values
// examples:
//   repeat_values([0, 1, 2], 1) -> [0, 1, 2]
//   repeat_values([0, 1, 2], 2) -> [0, 0, 1, 1, 2, 2]
//   repeat_values([0, 1, 2], 3) -> [0, 0, 0, 1, 1, 1, 2, 2, 2]
//   ...

template<typename Iterator>
class repeat_values
{
public:

	typedef typename thrust::iterator_difference<Iterator>::type difference_type;

	struct ediv_functor: public thrust::unary_function<difference_type,
			difference_type>
	{
		difference_type modulo;

		ediv_functor(difference_type modulo) :
			modulo(modulo)
		{
		}

		__host__ __device__
		difference_type operator()(const difference_type& i) const
		{
			return i / modulo;
		}
	};

	typedef typename thrust::counting_iterator<difference_type> CountingIterator;
	typedef typename thrust::transform_iterator<ediv_functor, CountingIterator> TransformIterator;
	typedef typename thrust::permutation_iterator<Iterator, TransformIterator> PermutationIterator;

	// type of the repeat_range iterator
	typedef PermutationIterator iterator;

	// construct repeat_range for the range [first,last)
	repeat_values(Iterator first, Iterator last, difference_type repeat) :
			first(first), last(last), repeat(repeat)
	{
	}

	iterator begin(void) const
	{
		return PermutationIterator(first,
				TransformIterator(CountingIterator(0), ediv_functor(last - first)));
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

#endif // __CUDA_REPEAT_VALUES_H__
