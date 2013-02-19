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
#include "CodeModulator_MFSK.h"
#include <assert.h>
#include <cstring>


//=================================================================================================
CodeModulator_MFSK::CodeModulator_MFSK(
		wsgc_float f_sampling,
		wsgc_float zero_frequency,
		wsgc_float frequency_shift,
		wsgc_float symbol_duration) :
	_local_oscillator(f_sampling, int(symbol_duration*f_sampling)),
	_zero_frequency(zero_frequency),
	_frequency_shift(frequency_shift),
	_symbol_duration(symbol_duration),
	_nb_symbol_samples(int(symbol_duration*f_sampling))
{}


//=================================================================================================
CodeModulator_MFSK::~CodeModulator_MFSK()
{}


//=================================================================================================
void CodeModulator_MFSK::modulate(wsgc_fftw_complex *out, std::vector<unsigned int>& symbols)
{
    std::vector<unsigned int>::const_iterator it = symbols.begin();
    const std::vector<unsigned int>::const_iterator it_end = symbols.end();

    for (unsigned int i=0; it != it_end; ++it, i++)
    {
    	_local_oscillator.make_next_samples(_zero_frequency + _frequency_shift*(*it));
    	memcpy(&out[_nb_symbol_samples*i], _local_oscillator.get_samples(), _nb_symbol_samples*sizeof(wsgc_complex));
    }
}

