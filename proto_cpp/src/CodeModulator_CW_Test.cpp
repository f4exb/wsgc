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
#include "CodeModulator_CW_Test.h"
#include <cstring>


void CodeModulator_CW_Test::fill_code_samples(wsgc_fftw_complex *fftw_code_in, std::vector<char>& code_bits)
{
    std::vector<char>::const_iterator it = code_bits.begin();
    const std::vector<char>::const_iterator it_end = code_bits.end();

    for (unsigned int i=0; it != it_end; ++it, i++)
    {
        // real part always on and imaginary part always off regardless of code
        fftw_code_in[i][0] = 1;
        fftw_code_in[i][1] = 0;    
    }
}


void CodeModulator_CW_Test::modulate(const wsgc_fftw_complex *in, wsgc_fftw_complex *out, std::vector<char>& code_bits)
{
    memcpy(out, in, code_bits.size()*sizeof(wsgc_fftw_complex));
}
