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

     MessageCorrelator

     Given the frequency and time displacement of correlation peak given by the
     Pilot Correlator it searches correlation for all PRNs of the message symbols
     to select the one that was sent. It uses straightforward time correlation.

*/

#include "MessageCorrelator_SoftPilot.h"
#include "GoldCodeGenerator.h"
#include "LocalCodes.h"
#include "WsgcUtils.h"
#include <cmath>
#include <cstring>

MessageCorrelator_SoftPilot::MessageCorrelator_SoftPilot(GoldCodeGenerator& gc_generator, LocalCodes& local_codes, wsgc_float f_sampling, wsgc_float f_chip, unsigned int prn_per_symbol, unsigned int nb_batch_prns) :
	_gc_generator(gc_generator),
	_local_codes(local_codes),
    _local_oscillator(f_sampling, gc_generator.get_nb_code_samples(f_sampling, f_chip)),
	_f_sampling(f_sampling),
	_f_chip(f_chip),
    _nb_msg_prns(_gc_generator.get_nb_message_codes()+1), // +1 for noise PRN
	_prn_per_symbol(prn_per_symbol),
	_fft_N(gc_generator.get_nb_code_samples(f_sampling,f_chip)),
    _delta_f(0.0),
    _correlation_matrices(_gc_generator.get_nb_message_codes()+1, prn_per_symbol, gc_generator.get_nb_code_samples(f_sampling,f_chip))
{
    static const wsgc_complex c_zero(0.0,0.0);

	// Allocate memory areas
	_src = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*sizeof(wsgc_fftw_complex));
	_corr_results = new wsgc_complex[_nb_msg_prns*_prn_per_symbol];
    _corr_delta_t = new unsigned int[_nb_msg_prns*_prn_per_symbol];
    _corr_avgsums = new wsgc_complex[_nb_msg_prns];

    for (unsigned int i=0; i<_nb_msg_prns; i++)
    {
        _corr_avgsums[i] = c_zero;
    
        for (unsigned int j=0; j<_prn_per_symbol; j++)
        {
            _corr_results[i*_prn_per_symbol + j] = c_zero;
        }
    }
}


MessageCorrelator_SoftPilot::~MessageCorrelator_SoftPilot()
{
	// Free memory areas
	 WSGC_FFTW_FREE(_src);
     delete[] _corr_avgsums;
     delete[] _corr_results;
}


void MessageCorrelator_SoftPilot::copy_input_samples(wsgc_complex *fft_src, wsgc_float delta_f)
{
    _delta_f = delta_f;

    if (delta_f == 0)
    {
        memcpy((void *) _src, (void *) fft_src, _fft_N*sizeof(wsgc_fftw_complex)); // no frequency correction, just copy
    }
    else
    {
        _local_oscillator.make_next_samples(delta_f);
        const wsgc_complex *lo_source_block = _local_oscillator.get_samples();
        
        for (unsigned int i=0; i<_fft_N; i++) // mix with LO to obtain zero IF 
        {
            _src[i] = fft_src[i] * lo_source_block[i];
        }
    }
}


void MessageCorrelator_SoftPilot::execute(unsigned int t_shift, unsigned int global_prn_index, CorrelationRecord& correlation_record)
{
    static const wsgc_complex c_zero(0.0,0.0);
    wsgc_complex corr_result;
    wsgc_float max_magnitude(0.0);
    wsgc_float avg_magnitude(0.0);
    wsgc_complex max_correlation;
    wsgc_float magnitude;
    wsgc_float max_noise;
    unsigned int max_prn_index(_gc_generator.get_nb_message_codes()); // init with noise PRN index
    const wsgc_complex *local_code;
    MessageCorrelationMatrices::CorrelationTuple_t correlation_tuple, noise_correlation_tuple;
    unsigned int pi = global_prn_index % _prn_per_symbol;
    
    for (unsigned int prni=0; prni < _nb_msg_prns; prni++) // search all message codes
    {
        local_code = _local_codes.get_local_code(prni);
        _corr_results[prni*_prn_per_symbol + pi] = c_zero;
        
        for (unsigned int i=0; i < _fft_N; i++) // do the time domain correlation
        {
            _corr_results[prni*_prn_per_symbol + pi] += _src[(t_shift + i) % _fft_N] * local_code[i];
        }
        
        _correlation_matrices.add_correlation_item(prni, pi, t_shift, _corr_results[prni*_prn_per_symbol + pi]);
        /* old process
        _corr_avgsums[prni] = _corr_results[prni*_prn_per_symbol];
        
        for (unsigned int ai=1; ai<_prn_per_symbol; ai++)
        {
            _corr_avgsums[prni] += _corr_results[prni*_prn_per_symbol + ai];
        }
        */
    }
        
    // new process
    _correlation_matrices.validate_prni_vector(pi);
    _correlation_matrices.process_averaging();
    _correlation_matrices.get_mag_max(correlation_tuple, max_magnitude, avg_magnitude);
    _correlation_matrices.get_noise_mag_max(noise_correlation_tuple, max_noise);

    // fill correlation record with correlation information part
    correlation_record.global_prn_index = global_prn_index;
    correlation_record.prn_per_symbol_index = global_prn_index % _prn_per_symbol;
    correlation_record.prn_index_max = correlation_tuple.prni;
    correlation_record.magnitude_max = max_magnitude;
    correlation_record.magnitude_avg = avg_magnitude;
    correlation_record.noise_max = max_noise;
    correlation_record.noise_avg = _correlation_matrices.get_noise_avg();
    correlation_record.shift_index_max = correlation_tuple.delta_t;
    correlation_record.pilot_shift = t_shift; // copy over
    correlation_record.f_rx = -_delta_f; // copy over
    correlation_record.phase_at_max = atan2(correlation_tuple.value.imag(), correlation_tuple.value.real());

    /* old process
    for (unsigned int prni=0; prni < _nb_msg_prns; prni++) // search all PRN averages to find maximum
    {
        WsgcUtils::magnitude_estimation(&_corr_avgsums[prni], &magnitude);
        
        if (prni == _nb_msg_prns - 1) // noise PRN
        {
        	max_noise = magnitude;
        }
        else
        {
			if (magnitude > max_magnitude)
			{
				max_magnitude = magnitude;
				max_correlation = _corr_avgsums[prni];
				max_prn_index = prni;
			}
        }
    }
    
    // fill correlation record with correlation information part
    correlation_record.global_prn_index = global_prn_index;
    correlation_record.prn_per_symbol_index = global_prn_index % _prn_per_symbol;
    correlation_record.prn_index_max = max_prn_index;
    correlation_record.magnitude_max = max_magnitude;
    correlation_record.noise_max = max_noise;
    correlation_record.shift_index_max = t_shift; // copy over
    correlation_record.f_rx = -_delta_f; // copy over
    correlation_record.phase_at_max = atan2(max_correlation.imag(), max_correlation.real());
    */
}

