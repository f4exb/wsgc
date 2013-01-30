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
      
     CodeModulator for DBPSK signals
     
*/
#include "CodeModulator_DBPSK.h"
#include "WsgcException.h"
#include <iostream>
#include <assert.h>


void CodeModulator_DBPSK::fill_code_samples(wsgc_fftw_complex *fftw_code_in, std::vector<char>& code_samples)
{
	throw WsgcException("fill_code_samples with code samples in input is not implemented for DBPSK");
}


void CodeModulator_DBPSK::modulate(const wsgc_fftw_complex *in, wsgc_fftw_complex *out, std::vector<char>& code_samples)
{
	throw WsgcException("modulate with code samples in input is not implemented for DBPSK");
}


void CodeModulator_DBPSK::fill_code_samples(wsgc_fftw_complex *fftw_code_in, std::vector<char>& code_bits, wsgc_float f_sampling, wsgc_float f_chip)
{
	wsgc_float fractional_index_increment;
	unsigned int number_of_samples = ((code_bits.size()/f_chip)*f_sampling);

    assert(!(f_sampling < f_chip));
    fractional_index_increment = ((wsgc_float) code_bits.size()) / number_of_samples;

    wsgc_float fractional_index = 0.0;
    int index;
    int index_1;
    char memory_bit = _seed_bit; // from last batch or init
    char code; // current code bit sample

    for (int i=0; i<number_of_samples; i++)
    {
        index = int(fractional_index);
        index_1 = int(fractional_index-fractional_index_increment);

        if ((index > index_1) || (i == 0)) // chip change
        {
        	if (code_bits[index] == 0) // phase change
        	{
        		code = (memory_bit ? 0 : 1);
        	}
        	else
        	{
        		code = memory_bit; // keep code bit
        	}
        }

        fftw_code_in[i][0] = 2*(code) - 1;
        fftw_code_in[i][1] = 2*(code) - 1;

        fractional_index += fractional_index_increment;
        memory_bit = code;
    }

    _seed_bit = code; // for next batch
}


void CodeModulator_DBPSK::modulate(const wsgc_fftw_complex *in, wsgc_fftw_complex *out, std::vector<char>& code_bits, wsgc_float f_sampling, wsgc_float f_chip)
{
	wsgc_float fractional_index_increment;
	unsigned int number_of_samples = ((code_bits.size()/f_chip)*f_sampling);

    assert(!(f_sampling < f_chip));
    fractional_index_increment = ((wsgc_float) code_bits.size()) / number_of_samples;

    wsgc_float fractional_index = 0.0;
    int index;
    int index_1;
    char memory_bit = _seed_bit; // from last batch or init
    char code; // current code bit sample

    for (int i=0; i<number_of_samples; i++)
    {
        index = int(fractional_index);
        index_1 = int(fractional_index-fractional_index_increment);

        if ((index > index_1) || (i == 0)) // chip change
        {
        	if (code_bits[index] == 0) // phase change
        	{
        		code = (memory_bit ? 0 : 1); // flip code bit
        	}
        	else
        	{
        		code = memory_bit; // keep code bit
        	}
        }

        out[i][0] = in[i][0] * (2*(code) - 1);
        out[i][1] = in[i][1] * (2*(code) - 1);

        fractional_index += fractional_index_increment;
        memory_bit = code;
    }

    _seed_bit = code; // for next batch
}
