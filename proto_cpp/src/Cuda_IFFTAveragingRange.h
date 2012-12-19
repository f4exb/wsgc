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

#ifndef __CUDA_IFFT_AVERAGING_RANGE_H__
#define __CUDA_IFFT_AVERAGING_RANGE_H__

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/functional.h>

/**
 * \brief Special iterator to make access for averaging to the strided values stored in the IFFT result
 *
 * This example illustrates how to make access for averaging to the strided values stored in the IFFT result
 *
 * Examples:
 *                                           even/odd (0/1) batch number.
 *                                 number of values in a batch (skip).  |
 *                              number of averaged values (repeat).  |  |
 *                                                                |  |  |
 *   ifft_averaging_range([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], 4, 3, 0) even batch:
 * [0, 1, 2, 3,  1, 2, 3, 4,  2, 3, 4, 5,  6, 7, 8, 9,  7, 8, 9, 10,  8, 9, 10, 11]
 *  0  1  2  3  4  5  6  7  8  9 10 11 
 *  +---------
 *     +---------
 *        +---------
 *                    +---------
 *                       +---------
 *                          +---------
 *
 *   ifft_averaging_range([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], 4, 3, 1) odd batch:
 * [3, 4, 5, 0,  4, 5, 0, 1,  5, 0, 1, 2,  9, 10, 11, 6,  10, 11, 6, 7,  11, 6, 7, 8]
 *  0  1  2  3  4  5  6  7  8  9 10 11 
 *  --       +-------
 *  ----        +----
 *  -------        +-
 *                    --       +-------
 *                    ----        +----
 *                    -------        +-
 */
template<typename Iterator>
class ifft_averaging_range
{
public:

	typedef typename thrust::iterator_difference<Iterator>::type difference_type;

	struct ifft_averaging_functor: public thrust::unary_function<difference_type,
			difference_type>
	{
		difference_type repeat;
		difference_type skip;
		difference_type odd;

		ifft_averaging_functor(difference_type repeat, difference_type skip, difference_type odd) :
				repeat(repeat), skip(skip), odd(odd)
		{
		}

		__host__ __device__
		difference_type operator()(const difference_type& i) const
		{
            return (((i%repeat) + ((i/repeat)%skip) + odd*skip)%(2*skip)) + (((i/repeat)/skip)*2*skip);
		}
	};

	typedef typename thrust::counting_iterator<difference_type> CountingIterator;
	typedef typename thrust::transform_iterator<ifft_averaging_functor, CountingIterator> TransformIterator;
	typedef typename thrust::permutation_iterator<Iterator, TransformIterator> PermutationIterator;

	// type of the ifft_averaging_range iterator
	typedef PermutationIterator iterator;

	// construct ifft_averaging_range for the range [first,last)
	ifft_averaging_range(Iterator first, Iterator last, difference_type repeat, difference_type skip, difference_type odd) :
			first(first), last(last), repeat(repeat), skip(skip), odd(odd), usable_length(((last-first)/(2*skip))*2*skip)
	{
	}

	iterator begin(void) const
	{
		return PermutationIterator(first,
				TransformIterator(CountingIterator(0), ifft_averaging_functor(repeat, skip, odd)));
	}

	iterator end(void) const
	{
		return begin() + ((usable_length/(2*skip))*skip*repeat); 
	}

protected:
	Iterator first;
	Iterator last;
	difference_type repeat;
	difference_type skip;
	difference_type odd;
	difference_type usable_length;
};


/**
 * \brief Special iterator to make access to the result of averaging the strided values stored in the IFFT result
 *
 * This example illustrates how to make access to the result of averaging the strided values stored in the IFFT result
 *
 * Examples:
 *                                          even/odd (0/1) batch number.
 *                                number of values in a batch (skip).  |
 *                             number of averaged values (repeat).  |  |
 *                                                               |  |  |
 *   ifft_averaged_range([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], 4, 3, 0) even batch:
 *   [0, 1, 2, 6, 7, 8]
 *
 *   ifft_averaged_range([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], 4, 3, 1) odd batch:
 *   [3, 4, 5, 9, 10, 11]
 */
template<typename Iterator>
class ifft_averaged_range
{
public:

	typedef typename thrust::iterator_difference<Iterator>::type difference_type;

	struct ifft_averaged_functor: public thrust::unary_function<difference_type,
			difference_type>
	{
		difference_type skip;
		difference_type odd;

		ifft_averaged_functor(difference_type skip, difference_type odd) :
				skip(skip), odd(odd)
		{
		}

		__host__ __device__
		difference_type operator()(const difference_type& i) const
		{
            return (i%skip) + odd*skip + (i/skip)*2*skip;
		}
	};

	typedef typename thrust::counting_iterator<difference_type> CountingIterator;
	typedef typename thrust::transform_iterator<ifft_averaged_functor, CountingIterator> TransformIterator;
	typedef typename thrust::permutation_iterator<Iterator, TransformIterator> PermutationIterator;

	// type of the ifft_averaged_range iterator
	typedef PermutationIterator iterator;

	// construct ifft_averaged_range for the range [first,last)
	ifft_averaged_range(Iterator first, Iterator last, difference_type skip, difference_type odd) :
			first(first), last(last), skip(skip), odd(odd), usable_length(((last-first)/(2*skip))*2*skip)
	{
	}

	iterator begin(void) const
	{
		return PermutationIterator(first,
				TransformIterator(CountingIterator(0), ifft_averaged_functor(skip, odd)));
	}

	iterator end(void) const
	{
		return begin() + ((usable_length/(2*skip))*skip); 
	}

protected:
	Iterator first;
	Iterator last;
	difference_type skip;
	difference_type odd;
	difference_type usable_length;
};

#endif // __CUDA_ifft_averaging_range_H__
