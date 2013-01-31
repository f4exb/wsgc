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
#include "WsgcUtils.h"
#include "CorrelationRecord.h"
#include "TrainingCorrelationRecord.h"
#include <iostream>
#include <assert.h>


const wsgc_complex DifferentialModulationMultiplePrnCorrelator_Host::_c_zero = (0.0, 0.0);

//=================================================================================================
DifferentialModulationMultiplePrnCorrelator_Host::DifferentialModulationMultiplePrnCorrelator_Host(
        wsgc_float f_sampling, 
        wsgc_float f_chip, 
		unsigned int prn_length,
        unsigned int prn_per_symbol,
        const std::vector<unsigned int>& prn_list,
        unsigned int symbol_window_size,
		std::vector<CorrelationRecord>& correlation_records,
		std::vector<TrainingCorrelationRecord>& training_correlation_records,
        const LocalCodesFFT_Host& local_codes_fft_base) :
        DifferentialModulationMultiplePrnCorrelator::DifferentialModulationMultiplePrnCorrelator(f_sampling, f_chip, prn_length, prn_per_symbol, prn_list, symbol_window_size, correlation_records, training_correlation_records),
        _local_codes_fft_base(local_codes_fft_base)
{
	// Input storage
    _samples = (wsgc_complex *) WSGC_FFTW_MALLOC((_int_samples_per_chip+(symbol_window_size*prn_per_symbol*_fft_N))*sizeof(wsgc_fftw_complex));
    _demod = (wsgc_complex *) WSGC_FFTW_MALLOC(symbol_window_size*prn_per_symbol*_fft_N*sizeof(wsgc_fftw_complex));
    _demod_fft = (wsgc_complex *) WSGC_FFTW_MALLOC(symbol_window_size*prn_per_symbol*_fft_N*sizeof(wsgc_fftw_complex));

    // IFFT items
    _ifft_in_tmp = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*sizeof(wsgc_fftw_complex));
    _ifft_out_tmp = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*sizeof(wsgc_fftw_complex));
    _corr_out = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*symbol_window_size*prn_per_symbol*prn_list.size()*sizeof(wsgc_fftw_complex));
    _corr_out_mag = new wsgc_float[_fft_N*symbol_window_size*prn_per_symbol*prn_list.size()];

    // Do the multiple FFT of input demodulated samples
    int N = _fft_N;
    int howmany = symbol_window_size*prn_per_symbol;
    _fft_plan = WSGC_FFTW_PLAN_MANY(1, &N, howmany,
    		reinterpret_cast<wsgc_fftw_complex *>(_demod),
    		0, 1, _fft_N,
    		reinterpret_cast<wsgc_fftw_complex *>(_demod_fft),
    		0, 1, _fft_N,
    		FFTW_FORWARD, FFTW_ESTIMATE);

    // Do the IFFT in temporary buffer
    _ifft_plan = WSGC_FFTW_PLAN(_fft_N,
                                reinterpret_cast<wsgc_fftw_complex *>(_ifft_in_tmp),
                                reinterpret_cast<wsgc_fftw_complex *>(_ifft_out_tmp),
                                FFTW_BACKWARD, FFTW_ESTIMATE);

    // Fill initial chip with zeros
    for (unsigned int i=0; i < _int_samples_per_chip; i++)
    {
    	_samples[i] = _c_zero;
    }

    _samples_length = _int_samples_per_chip;
}

        
//=================================================================================================
DifferentialModulationMultiplePrnCorrelator_Host::~DifferentialModulationMultiplePrnCorrelator_Host()
{
	// IFFT plan
	WSGC_FFTW_DESTROY_PLAN(_ifft_plan);

	// FFT plan
	WSGC_FFTW_DESTROY_PLAN(_fft_plan);

    // IFFT items
	delete[] _corr_out_mag;
	WSGC_FFTW_FREE(_corr_out);
	WSGC_FFTW_FREE(_ifft_out_tmp);
	WSGC_FFTW_FREE(_ifft_in_tmp);

	// Input storage
    WSGC_FFTW_FREE(_demod_fft);
    WSGC_FFTW_FREE(_demod);
    WSGC_FFTW_FREE(_samples);
}


//=================================================================================================
void DifferentialModulationMultiplePrnCorrelator_Host::set_initial_chip(wsgc_complex *chip_samples)
{
    memcpy(_samples, chip_samples, _int_samples_per_chip*sizeof(wsgc_complex));
    
    if (_samples_length == 0)
    {
        _samples_length = _int_samples_per_chip;
    }
}


//=================================================================================================
void DifferentialModulationMultiplePrnCorrelator_Host::set_samples(wsgc_complex *samples)
{
	assert(_samples_length+_fft_N <= _int_samples_per_chip+(_symbol_window_size*_prn_per_symbol*_fft_N));
    memcpy(&_samples[_samples_length], samples, _fft_N*sizeof(wsgc_complex));
    _samples_length += _fft_N;
    _global_prn_index++;
    _prns_length++;
}


//=================================================================================================
void DifferentialModulationMultiplePrnCorrelator_Host::execute_message()
{
	do_correlation();
	do_sum_averaging();
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
	differentially_demodulate_window();
	WSGC_FFTW_EXECUTE(_fft_plan);
	// for each PRN length in processing window half
	for (unsigned int prn_wi = 0; prn_wi < _symbol_window_size*_prn_per_symbol; prn_wi++)
	{
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
			_ifft_in_tmp[ffti] = _demod_fft[prn_wi*_fft_N + ffti] * _local_codes_fft_base.get_local_code(*prni_it)[ffti];
		}

		// do one IFFT
		WSGC_FFTW_EXECUTE(_ifft_plan);

		// push back the result in the global IFFT output array
		for (unsigned int ffti = 0; ffti < _fft_N; ffti++)
		{
			_corr_out[_symbol_window_size*_prn_per_symbol*_fft_N*prni_i + _symbol_window_size*_prn_per_symbol*ffti + prn_wi] = _ifft_out_tmp[ffti];
			//std::cout << prni_i << ":" << ffti << ":" << prn_wi << ":" << _symbol_window_size*_prn_per_symbol*_fft_N*prni_i + _symbol_window_size*_prn_per_symbol*ffti + prn_wi << std::endl;
		}
	}
}


//=================================================================================================
void DifferentialModulationMultiplePrnCorrelator_Host::do_sum_averaging()
{
	wsgc_float sumavg_mag, prni_max;

	init_results();

	for (unsigned int ci = 0; ci < _fft_N*_symbol_window_size*_prn_per_symbol*_prn_list.size(); ci++)
	{
		WsgcUtils::magnitude_estimation(&_corr_out[ci], &_corr_out_mag[ci]);
	}

	/*
	for (unsigned int prni_i = 0; prni_i < _prn_list.size(); prni_i++) // loop on correlated PRN numbers
	{
		prni_max = 0.0;

		for (unsigned int iffti = 0; iffti < _fft_N; iffti++) // for each time bin
		{
			for (unsigned int si = 0; si < _symbol_window_size; si++) // for each symbol
			{
				unsigned int start_symbol_iffti_index = _symbol_window_size*_prn_per_symbol*_fft_N*prni_i + _symbol_window_size*_prn_per_symbol*iffti + _prn_per_symbol*si;

				// calculate average complex sum on _prn_per_symbol samples
				for (unsigned int ai = 1; ai < _prn_per_symbol; ai++ )
				{
					_corr_out[start_symbol_iffti_index] += _corr_out[start_symbol_iffti_index + ai];
				}

				// take magnitude of (prni_i, si, iffti) result
				WsgcUtils::magnitude_estimation(&_corr_out[start_symbol_iffti_index], &sumavg_mag);

				// find biggest (prni_i, iffti) for the si
				if (sumavg_mag > _max_sy_mags[si])
				{
					_max_sy_mags[si] = sumavg_mag;
					_max_sy_iffti[si] = iffti;
					_max_sy_prni[si] = prni_i;
				}

				if (sumavg_mag > _max_sy_mags_prni[_symbol_window_size*prni_i+si])
				{
					_max_sy_mags_prni[_symbol_window_size*prni_i+si] = sumavg_mag;
				}
			}
		}
	}
	*/

	for (unsigned int prni_i = 0; prni_i < _prn_list.size(); prni_i++) // loop on correlated PRN numbers
	{
		prni_max = 0.0;

		for (unsigned int iffti = 0; iffti < _fft_N; iffti++) // for each time bin
		{
			for (unsigned int si = 0; si < _symbol_window_size; si++) // for each symbol
			{
				unsigned int start_symbol_iffti_index = _symbol_window_size*_prn_per_symbol*_fft_N*prni_i + _symbol_window_size*_prn_per_symbol*iffti + _prn_per_symbol*si;
				//std::cout << prni_i << ":" << iffti << ":" << si << ":" << start_symbol_iffti_index << std::endl;

				// calculate average sum on _prn_per_symbol samples
				for (unsigned int ai = 1; ai < _prn_per_symbol; ai++ )
				{
					_corr_out_mag[start_symbol_iffti_index] += _corr_out_mag[start_symbol_iffti_index + ai];
				}

				// find biggest (prni_i, iffti) for the si
				if (_corr_out_mag[start_symbol_iffti_index] > _max_sy_mags[si])
				{
					_max_sy_mags[si] = _corr_out_mag[start_symbol_iffti_index];
					_max_sy_iffti[si] = iffti;
					_max_sy_prni[si] = prni_i;
				}

				if (_corr_out_mag[start_symbol_iffti_index] > _max_sy_mags_prni[_symbol_window_size*prni_i+si])
				{
					_max_sy_mags_prni[_symbol_window_size*prni_i+si] = _corr_out_mag[start_symbol_iffti_index];
				}
			}
		}
	}

	// move results to correlation records
	for (unsigned int si = 0; si < _symbol_window_size; si++) // for each symbol
	{
		static const CorrelationRecord tmp_correlation_record(0);
		_correlation_records.push_back(tmp_correlation_record);
		CorrelationRecord& correlation_record = _correlation_records.back();

		correlation_record.global_prn_index = _global_prn_index - _symbol_window_size*_prn_per_symbol + (si+1)*_prn_per_symbol - 1; // point to last PRN of each symbol
		correlation_record.magnitude_max = _max_sy_mags[si];
		correlation_record.prn_index_max = _max_sy_prni[si];
		correlation_record.pilot_shift = _max_sy_iffti[si]; // there is no pilot but take the shift index place
		correlation_record.prn_per_symbol = _prn_per_symbol;
		correlation_record.prn_per_symbol_index = _prn_per_symbol-1;

		wsgc_float prn_max_sum = 0.0;
		wsgc_float noise_max;

		for (unsigned int prni_i = 0; prni_i < _prn_list.size(); prni_i++) // loop on correlated PRN numbers
		{
			if (prni_i < _prn_list.size()-1) // not noise PRN
			{
				prn_max_sum += _max_sy_mags_prni[_symbol_window_size*prni_i+si];
			}
			else // noise PRN
			{
				noise_max = _max_sy_mags_prni[_symbol_window_size*prni_i+si];
			}
		}

		correlation_record.magnitude_avg = (prn_max_sum - correlation_record.magnitude_max) / (_prn_list.size()-2);
		correlation_record.noise_max = noise_max;
	}
}


//=================================================================================================
void DifferentialModulationMultiplePrnCorrelator_Host::differentially_demodulate_window()
{
    for (unsigned int i=0; i < _symbol_window_size*_prn_per_symbol*_fft_N; i++)
    {
        _demod[i] = _samples[i] * std::conj(_samples[i+_int_samples_per_chip]);
    }
}


//=================================================================================================
void DifferentialModulationMultiplePrnCorrelator_Host::chip_carry()
{
    memcpy(_samples, &_samples[_symbol_window_size*_fft_N], _int_samples_per_chip*sizeof(wsgc_complex));
    _samples_length = _int_samples_per_chip;
    _prns_length = 0;
}
