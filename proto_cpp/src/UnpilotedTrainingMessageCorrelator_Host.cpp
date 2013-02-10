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

     This flavour of correlator deals with PRNs without use of a pilot sequence
     This is the Host version

*/
#include "UnpilotedTrainingMessageCorrelator_Host.h"
#include "WsgcUtils.h"
#include "TrainingCorrelationRecord.h"
#include <iostream>
#include <assert.h>


const wsgc_complex UnpilotedTrainingMessageCorrelator_Host::_c_zero = (0.0, 0.0);

//=================================================================================================
UnpilotedTrainingMessageCorrelator_Host::UnpilotedTrainingMessageCorrelator_Host(
        wsgc_float f_sampling, 
        wsgc_float f_chip, 
		unsigned int prn_length,
        unsigned int sequence_length,
        unsigned int averaging_length,
        const std::vector<unsigned int>& prn_list,
		std::vector<TrainingCorrelationRecord>& training_correlation_records,
        const LocalCodesFFT_Host& local_codes_fft_base) :
        UnpilotedTrainingMessageCorrelator::UnpilotedTrainingMessageCorrelator(f_sampling, f_chip, prn_length, sequence_length, averaging_length, prn_list, training_correlation_records),
        _local_codes_fft_base(local_codes_fft_base)
{
	assert(_sequence_length <= prn_list.size());

	// Input storage
    _samples = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*sizeof(wsgc_fftw_complex));
    _src_fft = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*sizeof(wsgc_fftw_complex));

    // IFFT items
    _ifft_in_tmp = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*sizeof(wsgc_fftw_complex));
    _ifft_out_tmp = (wsgc_complex *) WSGC_FFTW_MALLOC(_fft_N*sizeof(wsgc_fftw_complex));
	_corr_out_sums = new wsgc_complex[_fft_N*prn_list.size()];
    //_mag_avgsums = new wsgc_float[_fft_N*prn_list.size()];

    for (unsigned int i=0; i<prn_list.size()*_fft_N; i++)
    {
        //_mag_avgsums[i] = 0.0;
        _corr_out_sums[i] = _c_zero;
    }

    // Do the FFT of input demodulated samples
    _fft_plan = WSGC_FFTW_PLAN(_fft_N,
    		reinterpret_cast<wsgc_fftw_complex *>(_samples),
    		reinterpret_cast<wsgc_fftw_complex *>(_src_fft),
    		FFTW_FORWARD, FFTW_ESTIMATE);

    // Do the IFFT in temporary buffer
    _ifft_plan = WSGC_FFTW_PLAN(_fft_N,
                                reinterpret_cast<wsgc_fftw_complex *>(_ifft_in_tmp),
                                reinterpret_cast<wsgc_fftw_complex *>(_ifft_out_tmp),
                                FFTW_BACKWARD, FFTW_ESTIMATE);
}

        
//=================================================================================================
UnpilotedTrainingMessageCorrelator_Host::~UnpilotedTrainingMessageCorrelator_Host()
{
	// IFFT plan
	WSGC_FFTW_DESTROY_PLAN(_ifft_plan);

	// FFT plan
	WSGC_FFTW_DESTROY_PLAN(_fft_plan);

    // IFFT items
	//delete[] _mag_avgsums;
	delete[] _corr_out_sums;
	WSGC_FFTW_FREE(_ifft_out_tmp);
	WSGC_FFTW_FREE(_ifft_in_tmp);

	// Input storage
    WSGC_FFTW_FREE(_src_fft);
    WSGC_FFTW_FREE(_samples);
}


//=================================================================================================
void UnpilotedTrainingMessageCorrelator_Host::set_samples(wsgc_complex *samples)
{
	assert(_prn_in_seq_count < _sequence_length);

    memcpy(_samples, samples, _fft_N*sizeof(wsgc_complex));

    _prn_in_avg_count++;
    _prn_in_seq_count++;
    _global_prn_index++;
}


//=================================================================================================
void UnpilotedTrainingMessageCorrelator_Host::execute()
{
	do_correlation();

	if (_prn_in_avg_count == _averaging_length)
	{
		do_sum_averaging();
		_prn_in_avg_count = 0;

	    for (unsigned int i=0; i<_prn_list.size()*_fft_N; i++)
	    {
	        //_mag_avgsums[i] = 0.0;
	        _corr_out_sums[i] = _c_zero;
	    }
	}
}


//=================================================================================================
void UnpilotedTrainingMessageCorrelator_Host::do_correlation()
{
	unsigned int prni_i = 0;
	wsgc_float magnitude;

	WSGC_FFTW_EXECUTE(_fft_plan);

	for (unsigned int prni_i = _prn_in_avg_count-1; prni_i < _prn_list.size(); prni_i++) // loop on possible PRN numbers
	{
		// multiply source block by local code conjugate FFT
		for (unsigned int ffti = 0; ffti < _fft_N; ffti++) // multiply with local code
		{
			_ifft_in_tmp[ffti] = _src_fft[ffti] * _local_codes_fft_base.get_local_code(_prn_list[prni_i])[ffti];
		}

		// do one IFFT
		WSGC_FFTW_EXECUTE(_ifft_plan);

		// push back the result in the global IFFT output array
		for (unsigned int ffti = 0; ffti < _fft_N; ffti++)
		{
			//WsgcUtils::magnitude_estimation(&_ifft_out_tmp[ffti], &magnitude);
			_corr_out_sums[_fft_N*(prni_i-_prn_in_avg_count+1) + ffti] += _ifft_out_tmp[ffti];
			//_mag_avgsums[_fft_N*(prni_i-_prn_in_avg_count+1) + ffti] += magnitude;
		}
	}
}


//=================================================================================================
void UnpilotedTrainingMessageCorrelator_Host::do_sum_averaging()
{
	wsgc_float mag_max = 0.0;
	std::vector<wsgc_float> mag_sums;
	wsgc_float magnitude;
	unsigned int ci_max = 0;

	mag_sums.assign(_prn_list.size(), 0.0);

	for (unsigned int ci = 0; ci < _fft_N*_prn_list.size(); ci++)
	{
		WsgcUtils::magnitude_estimation(&_corr_out_sums[ci], &magnitude);
		mag_sums[ci/_fft_N] += magnitude;

		if (magnitude > mag_max)
		{
			mag_max = magnitude;
			ci_max = ci;
		}
	}

	static const TrainingCorrelationRecord tmp_training_correlation_record(_sequence_length, _averaging_length);
	_training_correlation_records.push_back(tmp_training_correlation_record);
	TrainingCorrelationRecord& correlation_record = _training_correlation_records.back();

	correlation_record._global_prn_index = _global_prn_index;
	correlation_record._prn_index_max = _prn_list[ci_max / _fft_N];
	correlation_record._pilot_shift = ci_max % _fft_N;
	correlation_record._magnitude_max = mag_max;
	correlation_record._magnitude_avgsum = (mag_sums[ci_max / _fft_N] - mag_max)/_prn_list.size();
}


