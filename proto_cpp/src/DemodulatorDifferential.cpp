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
      
     DemodulatorDifferential
      
*/
#include "DemodulatorDifferential.h"
#include "WsgcTypes.h"
#include <cstring>

//=================================================================================================
DemodulatorDifferential::DemodulatorDifferential(unsigned int differential_length) :
    _differential_length(differential_length)
{
    _differential_samples_at_origin = (wsgc_complex *) WSGC_FFTW_MALLOC(differential_length*sizeof(wsgc_fftw_complex));

    for (unsigned int i=0; i<_differential_length; i++)
    {
        _differential_samples_at_origin[i] = (0.0, 0.0);
    }
}


//=================================================================================================
DemodulatorDifferential::~DemodulatorDifferential()
{
    WSGC_FFTW_FREE(_differential_samples_at_origin);
}


//=================================================================================================
void DemodulatorDifferential::demodulate_in_place(wsgc_complex *samples_in, unsigned int samples_length)
{
    for (unsigned int i=0; i<_differential_length; i++)
    {
        samples_in[i] *= _differential_samples_at_origin[i];
    }

    for (unsigned int i=_differential_length; i < samples_length-_differential_length; i++)
    {
        samples_in[i] *= std::conj(samples_in[i+_differential_length]);
    }

    memcpy(_differential_samples_at_origin, &samples_in[samples_length-_differential_length], _differential_length*sizeof(wsgc_fftw_complex));
}


//=================================================================================================
void DemodulatorDifferential::demodulate_out_of_place(wsgc_complex *samples_in, wsgc_complex *samples_out, unsigned int samples_length)
{
    for (unsigned int i=0; i<_differential_length; i++)
    {
        samples_out[i] *= _differential_samples_at_origin[i];
    }

    for (unsigned int i=_differential_length; i < samples_length-_differential_length; i++)
    {
        samples_out[i] = samples_in[i] * std::conj(samples_in[i+_differential_length]);
    }
    
    memcpy(_differential_samples_at_origin, &samples_in[samples_length-_differential_length], _differential_length*sizeof(wsgc_fftw_complex));
}


//=================================================================================================
void DemodulatorDifferential::set_value_at_origin(wsgc_complex value)
{
    for (unsigned int i=0; i<_differential_length; i++)
    {
    	_differential_samples_at_origin[i] = value;
    }
}


//=================================================================================================
void DemodulatorDifferential::set_samples_at_origin(wsgc_complex *samples)
{
	memcpy(_differential_samples_at_origin, samples, _differential_length*sizeof(wsgc_fftw_complex));
}
