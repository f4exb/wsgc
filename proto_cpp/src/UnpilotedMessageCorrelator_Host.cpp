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

     UnpilotedMessageCorrelator

	 Does correlation on all possible message PRNs (normally also including noise PRN)
	 to find a correlation peak representing the time delay of the start of sequence.

	 Uses frequency domain correlation to check all possible delays at once which is
	 much more computationnaly efficient

	 This is the class for Host implementation

*/

#include "WsgcUtils.h"
#include "UnpilotedMessageCorrelator_Host.h"
#include "CorrelationRecord.h"
#include "CodeModulator.h"
#include "GoldCodeGenerator.h"
#include <string.h>


UnpilotedMessageCorrelator_Host::UnpilotedMessageCorrelator_Host(
		wsgc_float f_sampling,
		wsgc_float f_chip,
		unsigned int prn_per_symbol,
		unsigned int nb_batch_prns,
		std::vector<unsigned int>& message_symbols,
		CodeModulator& code_modulator,
		GoldCodeGenerator& gc_generator
		) :
	UnpilotedMessageCorrelator(f_sampling, f_chip, prn_per_symbol, nb_batch_prns),
	_local_codes(code_modulator, gc_generator, f_sampling, f_chip, message_symbols),
	_fft_N(gc_generator.get_nb_code_samples(f_sampling, f_chip)),
	_nb_msg_prns(_local_codes.get_nb_codes()),
	_prn_index(0),
	_batch_sum_magnitudes(nb_batch_prns*_nb_msg_prns,0.0)
{
	// FFT
    _fft_sample_in = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*sizeof(wsgc_complex));
    _fft_sample_out = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*sizeof(wsgc_complex));
    _fft_sample_plan = WSGC_FFTW_PLAN(_fft_N,
    		reinterpret_cast<wsgc_fftw_complex *>(_fft_sample_in),
    		reinterpret_cast<wsgc_fftw_complex *>(_fft_sample_out),
    		FFTW_FORWARD, FFTW_ESTIMATE);
    // IFFT
    _ifft_code_out = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*2*_nb_batch_prns*_nb_msg_prns*sizeof(wsgc_complex));
    _ifft_out_mags = (wsgc_float *) WSGC_FFTW_MALLOC(_fft_N*2*_nb_batch_prns*_nb_msg_prns*sizeof(wsgc_float));
    _ifft_code_in_tmp = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*sizeof(wsgc_complex));
    _ifft_code_out_tmp = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*sizeof(wsgc_complex));
    // Do the IFFT in temporary buffer
    _ifft_plan = WSGC_FFTW_PLAN(_fft_N,
    		reinterpret_cast<wsgc_fftw_complex *>(_ifft_code_in_tmp),
    		reinterpret_cast<wsgc_fftw_complex *>(_ifft_code_out_tmp),
    		FFTW_BACKWARD, FFTW_ESTIMATE);
}


UnpilotedMessageCorrelator_Host::~UnpilotedMessageCorrelator_Host()
{
    // IFFT plan
    WSGC_FFTW_DESTROY_PLAN(_ifft_plan);
    // IFFT items
    WSGC_FFTW_FREE(_ifft_code_out_tmp);
    WSGC_FFTW_FREE(_ifft_code_in_tmp);
    WSGC_FFTW_FREE(_ifft_code_out);
    // FFT plan
    WSGC_FFTW_DESTROY_PLAN(_fft_sample_plan);
    // FFT items
    WSGC_FFTW_FREE(_fft_sample_out);
    WSGC_FFTW_FREE(_fft_sample_in);
}


void UnpilotedMessageCorrelator_Host::set_source_block(wsgc_complex *source_block)
{
    memcpy(_fft_sample_in, source_block, _fft_N*sizeof(wsgc_complex));
}


bool UnpilotedMessageCorrelator_Host::execute(std::vector<CorrelationRecord>& message_correlation_records)
{
    bool averaging_done = false;
    wsgc_float mag;
    
    // do the FFT of the input signal

	WSGC_FFTW_EXECUTE(_fft_sample_plan);

    const wsgc_complex *local_code_block;

    for (unsigned int prni=0; prni < _nb_msg_prns; prni++)
    {
        // Multiply by the local code
    	local_code_block = _local_codes.get_local_code(prni);
        for (unsigned int ffti=0; ffti<_fft_N; ffti++)
        {
        	_ifft_code_in_tmp[ffti] = _fft_sample_out[ffti] * local_code_block[ffti];
        }

        // Do the IFFT in the temporary buffer
        WSGC_FFTW_EXECUTE(_ifft_plan);

        // Move result to batch location
        for (unsigned int iffti=0; iffti<_fft_N; iffti++)
        {
        	unsigned int batch_location = _prn_index % (2*_nb_batch_prns);
        	_ifft_code_out[2*_nb_batch_prns*_fft_N*prni + 2*_nb_batch_prns*iffti + batch_location] = _ifft_code_out_tmp[iffti];
        	WsgcUtils::magnitude_estimation(&_ifft_code_out_tmp[iffti], &mag);
        	_ifft_out_mags[2*_nb_batch_prns*_fft_N*prni + 2*_nb_batch_prns*iffti + batch_location] = mag;
        }
    }

    if ((_prn_index % _nb_batch_prns) == _nb_batch_prns-1) // last PRN in batch has been processed
    {
    	do_averaging(message_correlation_records);
        averaging_done = true;
    }

    _prn_index++;
    
    return averaging_done;
}


void UnpilotedMessageCorrelator_Host::do_averaging(std::vector<CorrelationRecord>& message_correlation_records)
{
    // _prn_index is at the last PRN processed, next falls on a 2*_nb_batch_prns boundary if this is an odd batch
	unsigned int batch_start = ((_prn_index+1)%(2*_nb_batch_prns) == 0 ? 0 : _nb_batch_prns); // this is slightly faster than below
	//unsigned int batch_start = (_prn_index+1)%(2*_nb_batch_prns); // should be same
	wsgc_complex z;
	wsgc_float sum_mag;
	_batch_max_magnitudes.assign(_nb_batch_prns, 0.0);
	_batch_sum_magnitudes.assign(_nb_batch_prns*_nb_msg_prns, 0.0);
	_batch_noise_max_magnitude.assign(_nb_batch_prns, 0.0);

    for (unsigned int prni=0; prni < _nb_msg_prns; prni++) // all PRNs
    {
    	for (unsigned int iffti=0; iffti<_fft_N; iffti++) // all IFFT positions
    	{
    		for (unsigned int bi=0; bi<_nb_batch_prns; bi++) // all positions in batch
    		{
    			/*
				for (unsigned int ai=1; ai<_prn_per_symbol; ai++) // sum average
				{
					_ifft_code_out[2*_nb_batch_prns*_fft_N*prni + 2*_nb_batch_prns*iffti + ((bi+batch_start)%(2*_nb_batch_prns))]
								   += _ifft_code_out[2*_nb_batch_prns*_fft_N*prni + 2*_nb_batch_prns*iffti + ((bi+batch_start+ai)%(2*_nb_batch_prns))];
				}
				*/

    			sum_mag = 0.0;

				for (unsigned int ai=0; ai<_prn_per_symbol; ai++) // sum average
				{
					sum_mag += _ifft_out_mags[2*_nb_batch_prns*_fft_N*prni + 2*_nb_batch_prns*iffti + ((bi+batch_start+ai)%(2*_nb_batch_prns))];
				}

				z = _ifft_code_out[2*_nb_batch_prns*_fft_N*prni + 2*_nb_batch_prns*iffti + ((bi+batch_start)%(2*_nb_batch_prns))]; // complex result of averaging at (prni, iffti, bi)
				//WsgcUtils::magnitude_estimation(&z, &sum_mag);
				//sum_mag = _ifft_out_mags[2*_nb_batch_prns*_fft_N*prni + 2*_nb_batch_prns*iffti + ((bi+batch_start)%(2*_nb_batch_prns))]; // average sum of magnitudes at (prni, iffti, bi)

				if (sum_mag > _batch_max_magnitudes[bi])
				{
					_batch_max_magnitudes[bi] = sum_mag;
					_batch_max_ifft_indexes[bi] = iffti;
					_batch_max_prn_indexes[bi] = prni;
					_batch_complex_values_max[bi] = z;
				}

				_batch_sum_magnitudes[bi*_nb_msg_prns+prni] += sum_mag;

				if (prni == _nb_msg_prns-1) // noise PRN
				{
					if (sum_mag > _batch_noise_max_magnitude[bi])
					{
						_batch_noise_max_magnitude[bi] = sum_mag;
					}
				}
    		}
    	}
    }

    static const CorrelationRecord tmp_correlation_record;

    for (unsigned int bi=0; bi<_nb_batch_prns; bi++)
    {
    	message_correlation_records.push_back(tmp_correlation_record);
    	CorrelationRecord& correlation_record = message_correlation_records.back();
    	correlation_record.global_prn_index = _prn_index - (_nb_batch_prns-1) + bi;
    	correlation_record.prn_per_symbol_index = correlation_record.global_prn_index % _prn_per_symbol;
    	correlation_record.prn_index_max = _batch_max_prn_indexes[bi];
    	correlation_record.magnitude_max = _batch_max_magnitudes[bi];
    	correlation_record.magnitude_avg = _batch_sum_magnitudes[bi*_nb_msg_prns+correlation_record.prn_index_max] / _fft_N;
    	correlation_record.shift_index_max = _batch_max_ifft_indexes[bi];
    	correlation_record.phase_at_max = atan2(_batch_complex_values_max[bi].imag(), _batch_complex_values_max[bi].real());
    	correlation_record.noise_max = _batch_noise_max_magnitude[bi];
    	correlation_record.noise_avg = _batch_sum_magnitudes[bi*_nb_msg_prns + _nb_msg_prns - 1] / _fft_N;
    	correlation_record.f_rx = 0.0; // N/A
    	correlation_record.frequency_locked = true; // N/A
    	correlation_record.pilot_shift = 0; // N/A
    	correlation_record.selected = true; // N/A
    }
}
    
