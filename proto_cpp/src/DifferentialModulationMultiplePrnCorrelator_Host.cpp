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

     DifferentialModulationMultiplePrnCorrelator

     This flavour of correlator deals with PRNs encoded with differential modulation
     This is the Host version

*/
#include "DifferentialModulationMultiplePrnCorrelator_Host.h"
#include <assert.h>

//=================================================================================================
DifferentialModulationMultiplePrnCorrelator_Host::DifferentialModulationMultiplePrnCorrelator_Host(
        wsgc_float f_sampling, 
        wsgc_float f_chip, 
		unsigned int prn_length,
        const std::vector<unsigned int>& prn_list,
        unsigned int prn_window_size,
        const LocalCodesFFT_Host& local_codes_fft_base) :
        DifferentialModulationMultiplePrnCorrelator::DifferentialModulationMultiplePrnCorrelator(f_sampling, f_chip, prn_length, prn_list, prn_window_size),
        _local_codes_fft_base(local_codes_fft_base)
{
	// Input storage
    _samples = (wsgc_complex *) WSGC_FFTW_MALLOC((_int_samples_per_chip+(prn_window_size*_fft_N))*sizeof(wsgc_fftw_complex));
    _demod = (wsgc_complex *) WSGC_FFTW_MALLOC(prn_window_size*_fft_N*sizeof(wsgc_fftw_complex));

    // IFFT items
    _ifft_in_tmp = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*sizeof(wsgc_fftw_complex));
    _ifft_out_tmp = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*sizeof(wsgc_fftw_complex));
    _corr_out = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*_prn_window_size*sizeof(wsgc_fftw_complex));

    // Do the IFFT in temporary buffer
    _ifft_plan = WSGC_FFTW_PLAN(_fft_N,
                                reinterpret_cast<wsgc_fftw_complex *>(_ifft_in_tmp),
                                reinterpret_cast<wsgc_fftw_complex *>(_ifft_out_tmp),
                                FFTW_BACKWARD, FFTW_ESTIMATE);

}

        
//=================================================================================================
DifferentialModulationMultiplePrnCorrelator_Host::~DifferentialModulationMultiplePrnCorrelator_Host()
{
	// IFFT plan
	WSGC_FFTW_DESTROY_PLAN(_ifft_plan);

    // IFFT items
	WSGC_FFTW_FREE(_corr_out);
	WSGC_FFTW_FREE(_ifft_out_tmp);
	WSGC_FFTW_FREE(_ifft_in_tmp);

	// Input storage
    WSGC_FFTW_FREE(_demod);
    WSGC_FFTW_FREE(_samples);
}


//=================================================================================================
void DifferentialModulationMultiplePrnCorrelator_Host::set_initial_chip(wsgc_complex *chip_samples)
{
    memcpy(_samples, chip_samples, _int_samples_per_chip);
    
    if (_samples_length == 0)
    {
        _samples_length = _int_samples_per_chip;
    }
}


//=================================================================================================
void DifferentialModulationMultiplePrnCorrelator_Host::set_samples(wsgc_complex *samples)
{
    assert(_samples_length+_fft_N < _int_samples_per_chip+(_prn_window_size*_fft_N));
    memcpy(&_samples[_samples_length], samples, _fft_N);
    _samples_length += _fft_N;
}


//=================================================================================================
void DifferentialModulationMultiplePrnCorrelator_Host::execute_message()
{
	do_correlation();
    // TODO: averaging
    chip_carry();
}


//=================================================================================================
void DifferentialModulationMultiplePrnCorrelator_Host::execute_training()
{
	do_correlation();
    // TODO: sliding averaging
    chip_carry();
}


//=================================================================================================
void DifferentialModulationMultiplePrnCorrelator_Host::do_correlation()
{
	for (unsigned int prn_wi = 0; prn_wi< _prn_window_size; prn_wi++) // for each PRN length in window
	{
		differentially_demodulate_window();
		do_correlation(prn_wi);
	}
}


//=================================================================================================
void DifferentialModulationMultiplePrnCorrelator_Host::do_correlation(unsigned int prn_wi)
{
	std::vector<unsigned int>::const_iterator prni_it = _prn_list.begin();
	const std::vector<unsigned int>::const_iterator prni_end = _prn_list.end();
	unsigned int prni_i = 0;

	for (; prni_it != prni_end; ++prni_it, prni_i++) // loop on possible PRN numbers
	{
		// multiply source block by local code conjugate FFT
		for (unsigned int ffti = 0; ffti < _fft_N; ffti++) // multiply with local code
		{
			_ifft_in_tmp[ffti] = _demod[prn_wi*_fft_N + ffti] * _local_codes_fft_base.get_local_code(*prni_it)[ffti];
		}

		// do one IFFT
		WSGC_FFTW_EXECUTE(_ifft_plan);

		// push back the result in the global IFFT output array
		for (unsigned int ffti = 0; ffti < _fft_N; ffti++)
		{
			_corr_out[_prn_window_size*ffti + prn_wi] = _ifft_out_tmp[ffti];
		}
	}
}


//=================================================================================================
void DifferentialModulationMultiplePrnCorrelator_Host::differentially_demodulate_window()
{
    for (unsigned int i=0; i < _prn_window_size*_fft_N; i++)
    {
        _demod[i] = _samples[i] * _samples[i+_int_samples_per_chip];
    }
}


//=================================================================================================
void DifferentialModulationMultiplePrnCorrelator_Host::chip_carry()
{
    memcpy(_samples, &_samples[_prn_window_size*_fft_N], _int_samples_per_chip);
}
