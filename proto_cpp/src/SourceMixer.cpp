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

     Source mixer

     Mixes two simulated source outputs
*/


#include "SourceMixer.h"
#include "SimulatedSource.h"
#include <algorithm>
#include <cmath>

SourceMixer::SourceMixer(SimulatedSource& source_A, SimulatedSource& source_B, wsgc_float b_gain) :
	_source_A(source_A),
	_source_B(source_B),
	_B_gain(b_gain),
	_nb_samples(std::max(_source_A.get_nb_samples(), _source_B.get_nb_samples()))
{
}

void SourceMixer::get_samples(wsgc_complex **samples)
{
    wsgc_complex *samples_A, *samples_B;
    wsgc_float normalization_factor = 1.0 / (1.0 + _B_gain);
	_source_A.create_samples(&samples_A);
	_source_B.create_samples(&samples_B);
	*samples = (wsgc_complex *) WSGC_FFTW_MALLOC(_nb_samples*sizeof(wsgc_fftw_complex));

	for (unsigned int i = 0; i < _nb_samples; i++)
	{
		if (i >= _source_A.get_nb_samples())
		{
		    (*samples)[i] = normalization_factor * _B_gain * samples_B[i];
		}
		else if (i >= _source_B.get_nb_samples())
		{
		    (*samples)[i] = normalization_factor * samples_A[i];
		}
		else
		{
		    (*samples)[i] = normalization_factor * (samples_A[i] + _B_gain * samples_B[i]);
		}
	}

	WSGC_FFTW_FREE(samples_A);
    WSGC_FFTW_FREE(samples_B);
}


SourceMixer::~SourceMixer()
{
}
