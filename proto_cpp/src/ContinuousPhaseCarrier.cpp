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
      
     Continuous phase carrier
      
     Makes pure sine carrier samples of specified length. The next batch of samples
     is in continuous phase with the previous
*/

#include <fftw3.h>
#include "ContinuousPhaseCarrier.h"

ContinuousPhaseCarrier::ContinuousPhaseCarrier(wsgc_float f_sampling, unsigned int length, wsgc_float init_phase) :
    _f_sampling(f_sampling),
    _length(length)
{
    _samples = (wsgc_complex *) WSGC_FFTW_MALLOC(_length*sizeof(wsgc_fftw_complex));
    _phase_mod_2pi = std::modf(init_phase, &_phase_intpart);
}


ContinuousPhaseCarrier::~ContinuousPhaseCarrier()
{
    WSGC_FFTW_FREE(_samples);
}

const wsgc_complex *ContinuousPhaseCarrier::get_samples() const
{
    return _samples;
}

void ContinuousPhaseCarrier::set_phase(wsgc_float phase)
{
    _phase_mod_2pi = std::modf(phase, &_phase_intpart);
}

void ContinuousPhaseCarrier::make_next_samples(wsgc_float f)
{
    wsgc_fftw_complex *fftw_samples = reinterpret_cast<wsgc_fftw_complex *>(_samples);
    
    for (int i=0; i<_length; i++)
    {
        fftw_samples[i][0] = cos(2.0 * M_PI * (((i*f)/_f_sampling) + _phase_mod_2pi)); // real
        fftw_samples[i][1] = sin(2.0 * M_PI * (((i*f)/_f_sampling) + _phase_mod_2pi)); // imaginary
    }    
    
    _phase_mod_2pi = std::modf(_phase_mod_2pi + ((_length*f)/_f_sampling), &_phase_intpart);
}
