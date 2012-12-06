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
      
     Source FFT - Host version
      
     Multiplies input samples by each of the sub-frequency LOs
     Do the FFT of the result
     
*/
#include "SourceFFT_Host.h"
#include "ContinuousPhaseCarrier.h"
#include <cstring>

SourceFFT_Host::SourceFFT_Host(wsgc_float f_sampling, 
                   wsgc_float f_chip,
                   unsigned int fft_N,
                   unsigned int freq_step_division) :
    SourceFFT(f_sampling, f_chip, fft_N, freq_step_division)
{
    _fft_sample_in  = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*_freq_step_division*sizeof(wsgc_fftw_complex));
    _fft_sample_out  = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*_freq_step_division*sizeof(wsgc_fftw_complex));
    
    for (unsigned int fi=0; fi < _freq_step_division; fi++)
    {
        _fft_sample_plans.push_back(WSGC_FFTW_PLAN(_fft_N,
                                                   reinterpret_cast<wsgc_fftw_complex *>(_fft_sample_in + fi*_fft_N),
                                                   reinterpret_cast<wsgc_fftw_complex *>(_fft_sample_out + fi*_fft_N),
                                                   FFTW_FORWARD, FFTW_ESTIMATE)
                                   );
    	_local_oscillators.push_back(new ContinuousPhaseCarrier(_f_sampling, _fft_N));
    }
}


SourceFFT_Host::~SourceFFT_Host()
{
    // Local oscillators
	for (std::vector<ContinuousPhaseCarrier*>::iterator it=_local_oscillators.begin(); it != _local_oscillators.end(); ++it)
	{
		delete *it;
	}

    // FFT Plans
    std::vector<wsgc_fftw_plan>::iterator fft_plan_it = _fft_sample_plans.begin();
    std::vector<wsgc_fftw_plan>::const_iterator fft_plan_end = _fft_sample_plans.end();

    for (; fft_plan_it != fft_plan_end; ++fft_plan_it)
    {
        WSGC_FFTW_DESTROY_PLAN(*fft_plan_it);
    }
    
    // FFT items
    WSGC_FFTW_FREE(_fft_sample_in);
    WSGC_FFTW_FREE(_fft_sample_out);
}


const wsgc_complex *SourceFFT_Host::get_fft_samples(wsgc_complex *source_block) const
{
	wsgc_float freq_interstep = _f_sampling / (_fft_N * _freq_step_division);

    for (unsigned int fsi=0; fsi < _freq_step_division; fsi++)
    {
    	if (fsi == 0)
    	{
    		memcpy((void *) _fft_sample_in, (const void *) source_block, _fft_N*sizeof(wsgc_fftw_complex)); // at f=0 just copy input samples to input buffer
    	}
    	else
    	{
    		_local_oscillators[fsi-1]->make_next_samples(fsi*freq_interstep);
			 const wsgc_complex *lo_source_block = _local_oscillators[fsi-1]->get_samples();

    		for(unsigned int i=0; i<_fft_N; i++)
    		{
    			_fft_sample_in[fsi*_fft_N + i] = source_block[i] * lo_source_block[i];
    		}
    	}

    	WSGC_FFTW_EXECUTE(_fft_sample_plans[fsi]);
    }
    
    return _fft_sample_out;
}
