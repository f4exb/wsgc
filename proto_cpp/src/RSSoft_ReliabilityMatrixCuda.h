/*
     Copyright 2012-2013 Edouard Griffiths <f4exb at free dot fr>

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

     RSSoft_ReliabilityMatrixCuda.h

     Cuda implementation of a reliability matrix in device memory
*/
#ifndef __RSSOFT_RELIABILITY_MATRIX_CUDA_H__
#define __RSSOFT_RELIABILITY_MATRIX_CUDA_H__

#include <thrust/device_vector.h>

namespace rssoft
{
	class ReliabilityMatrix;
}


class RSSoft_ReliabilityMatrixCuda
{
public:
	/**
	 * Constructor
	 * \param nb_symbols_log2 Log2 of the number of symbols used (number of symbols is a power of two)
	 * \param message_length Length of one message block to be decoded
	 */
	RSSoft_ReliabilityMatrixCuda(unsigned int nb_symbols_log2, unsigned int message_length);

	/**
	 * Destructor. Frees the matrix storage.
	 */
	~RSSoft_ReliabilityMatrixCuda();

	/**
	 * Enter one more symbol position data
	 * \param symbol_data Pointer to symbol data array. There must be nb_symbol values corresponding to the relative reliability of each symbol for the current symbol position in the message
	 */
	void enter_symbol_data(thrust::device_vector<float> d_symbol_data);

    /**
     * Enter an erasure at current symbol position. This is done by zeroing out the corresponding column in the matrix thus neutralizing it for further multiplicity calculation.
     */
    void enter_erasure();

    /**
     * Copy to a RSSoft library's reliability matrix internal storage on host
     */
    void copy_to_host(rssoft::ReliabilityMatrix& reliability_matrix);

	/**
	 * Resets the message symbol counter
	 */
	void reset_message_symbol_count()
	{
		message_symbol_count = 0;
	}

	/**
	 * Get the log2 of the number of symbols (i.e. rows)
	 */
	unsigned int get_nb_symbols_log2() const
	{
		return nb_symbols_log2;
	}

	/**
	 * Get the number of symbols (i.e. rows)
	 */
	unsigned int get_nb_symbols() const
	{
		return nb_symbols;
	}

	/**
	 * Get the number of message symbols (i.e. columns)
	 */
	unsigned int get_message_length() const
	{
		return message_length;
	}

protected:
	unsigned int nb_symbols_log2;
	unsigned int nb_symbols;
	unsigned int message_length;
	unsigned int message_symbol_count; //!< incremented each time a new message symbol data is entered
	thrust::device_vector<float> matrix; //!< The reliability matrix stored column first
};

#endif // __RSSOFT_RELIABILITY_MATRIX_CUDA_H__
