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

     Source FFT - CUDA version

     Multiplies input samples by each of the sub-frequency LOs
     Do the FFT of the result

*/
#include "SourceFFT_Cuda.h"
#include "ContinuousPhaseCarrier.h"
#include "WsgcException.h"

SourceFFT_Cuda::SourceFFT_Cuda(wsgc_float f_sampling,
                   wsgc_float f_chip,
                   unsigned int fft_N,
                   unsigned int freq_step_division) :
    SourceFFT(f_sampling, f_chip, fft_N, freq_step_division),
    _d_fft_in(_fft_N*_freq_step_division),
    _d_fft_out(_fft_N*_freq_step_division)
{

    _h_fft_sample_in  = new wsgc_complex[_fft_N*_freq_step_division];
    cufftResult_t fft_stat = cufftPlan1d(&_fft_plan, _fft_N, CUFFT_C2C, _freq_step_division);
    //TODO: check return code

    for (unsigned int fi=0; fi < _freq_step_division; fi++)
    {
    	_local_oscillators.push_back(new ContinuousPhaseCarrier(_f_sampling, _fft_N));
    }
}


SourceFFT_Cuda::~SourceFFT_Cuda()
{
    // Local oscillators
	for (std::vector<ContinuousPhaseCarrier*>::iterator it=_local_oscillators.begin(); it != _local_oscillators.end(); ++it)
	{
		delete *it;
	}

    cufftDestroy(_fft_plan); // FFT Plan
    delete[] _h_fft_sample_in; // Host FFT items
}


const thrust::device_vector<cuComplex>& SourceFFT_Cuda::do_fft(wsgc_complex *source_block)
{
    wsgc_float freq_interstep = _f_sampling / (_fft_N * _freq_step_division);

    // multiply by sub-frequency LOs

    for (unsigned int fsi=0; fsi < _freq_step_division; fsi++)
    {
    	if (fsi == 0)
    	{
    		memcpy((void *) _h_fft_sample_in, (const void *) source_block, _fft_N*sizeof(wsgc_fftw_complex)); // at f=0 just copy input samples to input buffer
    	}
    	else
    	{
    		_local_oscillators[fsi-1]->make_next_samples(fsi*freq_interstep);
			 const wsgc_complex *lo_source_block = _local_oscillators[fsi-1]->get_samples();

    		for(unsigned int i=0; i<_fft_N; i++)
    		{
    			_h_fft_sample_in[fsi*_fft_N + i] = source_block[i] * lo_source_block[i];
    		}
    	}
    }

    // copy to device vector

    thrust::copy(
        reinterpret_cast<const cuComplex *>(_h_fft_sample_in),
        reinterpret_cast<const cuComplex *>(_h_fft_sample_in+_fft_N*_freq_step_division),
        _d_fft_in.begin()
    );

    // Do FFT

    if (cufftExecC2C(_fft_plan, thrust::raw_pointer_cast(&_d_fft_in[0]), thrust::raw_pointer_cast(&_d_fft_out[0]), CUFFT_FORWARD) != CUFFT_SUCCESS)
    {
        throw WsgcException("CUFFT Error: Failed to do FFT of source");
    }

    return _d_fft_out;
}
