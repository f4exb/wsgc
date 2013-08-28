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
#include "RSSoft_ReliabilityMatrixCuda.h"
#include "RS_ReliabilityMatrix.h"

//=================================================================================================
RSSoft_ReliabilityMatrixCuda::RSSoft_ReliabilityMatrixCuda(unsigned int _nb_symbols_log2, unsigned int _message_length) :
	nb_symbols_log2(_nb_symbols_log2),
	nb_symbols(1<<_nb_symbols_log2),
	message_length(_message_length),
	message_symbol_count(0),
	matrix(nb_symbols*message_length)
{}

//=================================================================================================
RSSoft_ReliabilityMatrixCuda::~RSSoft_ReliabilityMatrixCuda()
{}

//=================================================================================================
void RSSoft_ReliabilityMatrixCuda::enter_symbol_data(thrust::device_vector<float> d_symbol_data)
{
	if (message_symbol_count < message_length)
	{
        thrust::copy(d_symbol_data.begin(), d_symbol_data.end(), matrix.begin() + message_symbol_count*nb_symbols);
		message_symbol_count++;
	}
}

//=================================================================================================
void RSSoft_ReliabilityMatrixCuda::enter_erasure()
{
	if (message_symbol_count < message_length)
	{
        thrust::fill(matrix.begin() + message_symbol_count*nb_symbols, matrix.begin() + (message_symbol_count+1)*nb_symbols, 0);
    }
}

//=================================================================================================
void RSSoft_ReliabilityMatrixCuda::copy_to_host(rssoft::RS_ReliabilityMatrix& reliability_matrix)
{
	float *host_matrix = (float *) reliability_matrix.get_raw_matrix();
    thrust::copy(matrix.begin(), matrix.end(), host_matrix);
    message_symbol_count = 0;
}

