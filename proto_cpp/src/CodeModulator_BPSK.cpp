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
      
     CodeModulator
     
*/
#include "CodeModulator_BPSK.h"
#include <assert.h>


void CodeModulator_BPSK::fill_code_samples(wsgc_fftw_complex *fftw_code_in, std::vector<char>& code_bits)
{
    std::vector<char>::const_iterator it = code_bits.begin();
    const std::vector<char>::const_iterator it_end = code_bits.end();

    for (unsigned int i=0; it != it_end; ++it, i++)
    {
        // +/- 180 degrees with unit module (BPSK)
        fftw_code_in[i][0] = 2*(*it) - 1; // 0/1 -> -1/+1 
        fftw_code_in[i][1] = 0;    
    }
}


void CodeModulator_BPSK::fill_code_samples(wsgc_fftw_complex *fftw_code_in, std::vector<char>& code_bits, wsgc_float f_sampling, wsgc_float f_chip)
{
	wsgc_float fractional_index_increment;
	unsigned int number_of_samples = ((code_bits.size()/f_chip)*f_sampling);

    assert(!(f_sampling < f_chip));
    fractional_index_increment = ((wsgc_float) code_bits.size()) / number_of_samples;

    wsgc_float fractional_index = 0.0;
    unsigned int index = 0;

    for (int i=0; i<number_of_samples; i++)
    {
        index = int(fractional_index);
        fftw_code_in[i][0] = 2*(code_bits[index]) - 1;
        fftw_code_in[i][1] = 0;

        fractional_index += fractional_index_increment;
    }
}


void CodeModulator_BPSK::modulate(const wsgc_fftw_complex *in, wsgc_fftw_complex *out, std::vector<char>& code_bits)
{
    std::vector<char>::const_iterator it = code_bits.begin();
    const std::vector<char>::const_iterator it_end = code_bits.end();

    for (unsigned int i=0; it != it_end; ++it, i++)
    {
        out[i][0] = in[i][0] * (2*(*it)-1);
        out[i][1] = in[i][1] * (2*(*it)-1);
    }
}


void CodeModulator_BPSK::modulate(const wsgc_fftw_complex *in, wsgc_fftw_complex *out, std::vector<char>& code_bits, wsgc_float f_sampling, wsgc_float f_chip)
{
	wsgc_float fractional_index_increment;
	unsigned int number_of_samples = ((code_bits.size()/f_chip)*f_sampling);

    assert(!(f_sampling < f_chip));
    fractional_index_increment = ((wsgc_float) code_bits.size()) / number_of_samples;

    wsgc_float fractional_index = 0.0;
    unsigned int index = 0;

    for (int i=0; i<number_of_samples; i++)
    {
        index = int(fractional_index);
        out[i][0] = in[i][0] * (2*(code_bits[index]) - 1);
        out[i][1] = in[i][1] * (2*(code_bits[index]) - 1);

        fractional_index += fractional_index_increment;
    }
}
