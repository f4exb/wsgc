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

#include "PilotedMessageCorrelator_Cuda.h"
#include "GoldCodeGenerator.h"
#include "PilotCorrelationAnalyzer.h"
#include "LocalCodes_Cuda.h"
#include "WsgcUtils.h"
#include "Cuda_Operators.h"
#include "Cuda_RepeatRange.h"
#include "Cuda_RepeatValue.h"
#include "Cuda_ShiftedBySegmentsRange.h"
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/transform_reduce.h>
#include <cmath>
#include <cstring>

PilotedMessageCorrelator_Cuda::PilotedMessageCorrelator_Cuda(
		LocalCodes_Cuda& local_codes,
		wsgc_float f_sampling,
		wsgc_float f_chip,
		unsigned int prn_per_symbol) :
	PilotedMessageCorrelator(f_sampling, f_chip, prn_per_symbol),
	_local_codes(local_codes),
    _local_oscillator(f_sampling, _local_codes.get_gc_generator().get_nb_code_samples(f_sampling, f_chip)),
    _nb_msg_prns(_local_codes.get_gc_generator().get_nb_message_codes()+1), // +1 for noise PRN
	_fft_N(_local_codes.get_gc_generator().get_nb_code_samples(f_sampling,f_chip)),
	_d_corr_in(_fft_N),
	_d_mul_out(_fft_N*_nb_msg_prns),
    _d_keys(_nb_msg_prns),
    _d_corr_out(_nb_msg_prns)
{
	// Allocate memory areas
    _src = new wsgc_complex[_fft_N];
}


PilotedMessageCorrelator_Cuda::~PilotedMessageCorrelator_Cuda()
{
	// Free memory areas
    delete[] _src;
}


void PilotedMessageCorrelator_Cuda::execute(PilotCorrelationAnalyzer& pilot_correlation_analyzer)
{
    float avgsum_magnitude;
    unsigned int pi = pilot_correlation_analyzer.get_start_analysis_pilot_correlation_records_index();
    const std::vector<PilotCorrelationRecord> pilot_correlation_records = pilot_correlation_analyzer.get_pilot_correlation_records();
    int cuda_nb_msg_prns = _nb_msg_prns;

    for (unsigned int pai=0; pai < pilot_correlation_analyzer.get_analysis_window_size()*_prn_per_symbol; pi++, pai++) // scan PRNs in the window
    {
    	CorrelationRecord& correlation_record = pilot_correlation_analyzer.new_message_correlation_record(pilot_correlation_records[pi].prn_index);
    	float delta_f = pilot_correlation_records[pi].delta_f;
    	unsigned int delta_t = pilot_correlation_records[pi].t_index_max;
    	std::complex<float> *prn_src = pilot_correlation_analyzer.get_samples(pilot_correlation_records[pi].prn_index);

    	// make input samples if necessary
    	if (pilot_correlation_records[pi].selected)
    	{
        	// make input samples by multiplying PRN samples by local oscillator samples
    		if (delta_f == 0.0)
    		{
    			memcpy((void *) _src, (void *) prn_src, _fft_N*sizeof(wsgc_fftw_complex)); // no frequency correction, just copy
    		}
    		else
    		{
    			 _local_oscillator.make_next_samples(delta_f);
    			 const wsgc_complex *lo_source_block = _local_oscillator.get_samples();

				for (unsigned int i=0; i < _fft_N; i++)
				{
					_src[i] = prn_src[i] * lo_source_block[i];
				}
    		}

    		correlation_record.selected = true;
    	}
    	else
    	{
    		correlation_record.selected = false;
    	}

    	// copy multiplied source to input device vector
    	//_d_corr_in.assign(reinterpret_cast<const cuComplex *>(_src), reinterpret_cast<const cuComplex *>(_src+_fft_N));
        thrust::copy(
            reinterpret_cast<const cuComplex *>(_src),
            reinterpret_cast<const cuComplex *>(_src+_fft_N),
            _d_corr_in.begin()
        );

        // create a repetition of the source for all PRNs
        repeat_range<thrust::device_vector<cuComplex>::iterator> d_corr_in_multi(_d_corr_in.begin(), _d_corr_in.end(), _nb_msg_prns);

        // multiply by the local codes

        const thrust::device_vector<cuComplex>& d_codes = _local_codes.get_local_codes(); // get the serial matrix of local codes
        shifted_by_segments_range<thrust::device_vector<cuComplex>::const_iterator> d_codes_shifted(d_codes.begin(), d_codes.end(), _fft_N, -delta_t); // shift the local codes
        thrust::transform(d_codes_shifted.begin(), d_codes_shifted.end(), d_corr_in_multi.begin(), _d_mul_out.begin(), cmulc_functor2()); // multiply together with many copies of source

        // reduce by key to get correlation value (sum of _d_mul_out segments) for each PRN, result keys will be PRN numbers incidentally
        
        repeat_values<thrust::counting_iterator<int> > key_counter(thrust::make_counting_iterator(0), thrust::make_counting_iterator(cuda_nb_msg_prns), _fft_N);
        thrust::reduce_by_key(key_counter.begin(), key_counter.end(), _d_mul_out.begin(), _d_keys.begin(), _d_corr_out.begin(), thrust::equal_to<int>(), caddc_functor());
        
        // TODO: missing averaging on prn_per_symbol. This should be a prefix sum + stride kind of pb
        
        // further reduce the result (values) to get the sum of all magnitudes in order to compute average
        avgsum_magnitude = thrust::transform_reduce(_d_corr_out.begin(), _d_corr_out.end(), mag_squared_functor<cuComplex, float>(), 0.0, thrust::plus<float>());
        // search PRN index with largest magnitude (squared norm)
        thrust::device_vector<cuComplex>::iterator max_prn_it = thrust::max_element(_d_corr_out.begin(), _d_corr_out.end(), lesser_mag_squared<cuComplex>());    
        
        cuComplex z = *max_prn_it;

    	correlation_record.prn_per_symbol_index = pai % _prn_per_symbol;
    	correlation_record.prn_index_max = max_prn_it - _d_corr_out.begin();
    	correlation_record.magnitude_max = mag_squared_functor<cuComplex, float>()(z);
    	correlation_record.magnitude_avg = avgsum_magnitude / (_nb_msg_prns-1);
    	correlation_record.noise_avg = mag_squared_functor<cuComplex, float>()(_d_corr_out.back()); // TODO: as averaging is not implemented yet this is not actually an average
    	correlation_record.pilot_shift = delta_t;
    	correlation_record.shift_index_max = delta_t; // copy over
    	correlation_record.f_rx = -delta_f; // copy over
    	correlation_record.phase_at_max = atan2(z.y, z.x);
    }
}

