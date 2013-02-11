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

     This is the CUDA implementation
     
*/

#include "PilotedMessageCorrelator_Cuda.h"
#include "GoldCodeGenerator.h"
#include "PilotCorrelationAnalyzer.h"
#include "LocalCodes_Cuda.h"
#include "WsgcUtils.h"
#include "Cuda_Operators.h"
#include "Cuda_RepeatRange.h"
#include "Cuda_RepeatValue.h"
#include "Cuda_StridedShiftedRange.h"
#include "Cuda_StridedFoldedRange.h"
#include "Cuda_StridedRange.h"
#include "Cuda_ShiftedBySegmentsRange.h"
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/transform_reduce.h>
#include <thrust/scan.h>
#include <cmath>
#include <cstring>

const cuComplex PilotedMessageCorrelator_Cuda::_c_zero = {0.0, 0.0};
const PilotedMessageCorrelator_Cuda::transient_corr_value_t PilotedMessageCorrelator_Cuda::_init_transient_corr_value = {false, 0, 0.0, 0.0, 0.0};

PilotedMessageCorrelator_Cuda::PilotedMessageCorrelator_Cuda(
		LocalCodes_Cuda& local_codes,
		wsgc_float f_sampling,
		wsgc_float f_chip,
		unsigned int prn_per_symbol,
		unsigned int cuda_device) :
	CudaDeviceManager::CudaDeviceManager(cuda_device),
	PilotedMessageCorrelator(f_sampling, f_chip, prn_per_symbol),
	_local_codes(local_codes),
    _local_oscillator(f_sampling, _local_codes.get_gc_generator().get_nb_code_samples(f_sampling, f_chip)),
    _nb_msg_prns(_local_codes.get_gc_generator().get_nb_message_codes()+1), // +1 for noise PRN
	_fft_N(_local_codes.get_gc_generator().get_nb_code_samples(f_sampling,f_chip)),
	_d_corr_in(_fft_N),
	_d_mul_out(_fft_N*_nb_msg_prns),
    _d_corr_out(_nb_msg_prns*_prn_per_symbol),
    _d_corr_mag(_nb_msg_prns*_prn_per_symbol),
    _d_corr_mag_avgsum(_nb_msg_prns*_prn_per_symbol),
    _h_corr_mag_avgsum(_nb_msg_prns*_prn_per_symbol),
	_d_keys(_nb_msg_prns)
{
	thrust::fill(_d_corr_out.begin(), _d_corr_out.end(), _c_zero);
	// Allocate memory areas
    _src = new wsgc_complex[_fft_N];
    _transient_corr_values.assign(_prn_per_symbol, _init_transient_corr_value);
    _max_avg.assign(_prn_per_symbol, 0.0);
    _max_avg_index.assign(_prn_per_symbol, 0);
}


PilotedMessageCorrelator_Cuda::~PilotedMessageCorrelator_Cuda()
{
	// Free memory areas
    delete[] _src;
}


void PilotedMessageCorrelator_Cuda::execute(PilotCorrelationAnalyzer& pilot_correlation_analyzer)
{
    //float avgsum_magnitude;
    unsigned int pi = pilot_correlation_analyzer.get_start_analysis_pilot_correlation_records_index();
    const std::vector<PilotCorrelationRecord> pilot_correlation_records = pilot_correlation_analyzer.get_pilot_correlation_records();
    //int cuda_nb_msg_prns = _nb_msg_prns;

    for (unsigned int pai=0; pai < pilot_correlation_analyzer.get_analysis_window_size()*_prn_per_symbol; pi++, pai++) // scan PRNs in the window
    {
    	//CorrelationRecord& correlation_record = pilot_correlation_analyzer.new_message_correlation_record(pilot_correlation_records[pi].prn_index);
    	float delta_f = pilot_correlation_records[pi].delta_f;
    	unsigned int delta_t = pilot_correlation_records[pi].t_index_max;

        _transient_corr_values[pai % _prn_per_symbol].selected = pilot_correlation_records[pi].selected;
        _transient_corr_values[pai % _prn_per_symbol].global_prn_i = pilot_correlation_records[pi].prn_index;
        _transient_corr_values[pai % _prn_per_symbol].delta_t = delta_t;
        _transient_corr_values[pai % _prn_per_symbol].delta_f = delta_f;
        _transient_corr_values[pai % _prn_per_symbol].pilot_phase_at_max = pilot_correlation_records[pi].phase_at_max;

    	std::complex<float> *prn_src = pilot_correlation_analyzer.get_samples(pilot_correlation_records[pi].prn_index);

        // zero out _d_corr_out array at each start of symbol
        if (pai % _prn_per_symbol == 0)
        {
            thrust::fill(_d_corr_out.begin(), _d_corr_out.end(), _c_zero);
        }
        
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

			repeat_values<thrust::counting_iterator<int> > key_counter(thrust::make_counting_iterator(0), thrust::make_counting_iterator((int)_nb_msg_prns), _fft_N);
			strided_shifted_range<thrust::device_vector<cuComplex>::iterator> d_corr_out_avg_in(_d_corr_out.begin(), _d_corr_out.end(), _prn_per_symbol, pai%_prn_per_symbol);
			thrust::reduce_by_key(key_counter.begin(), key_counter.end(), _d_mul_out.begin(), _d_keys.begin(), d_corr_out_avg_in.begin(), thrust::equal_to<int>(), caddc_functor());
    	}

    	// calculate averaging sums etc... for whole symbol at once

    	if (pai % _prn_per_symbol == _prn_per_symbol-1)
    	{
            // calculate magnitudes (squared norm)
            thrust::transform(_d_corr_out.begin(), _d_corr_out.end(), _d_corr_mag.begin(), mag_squared_functor<cuComplex, float>());
            
			// progressive averaging sum (that is inclusive scan) of magnitudes for each PRN. This gives the successive averaging sums for each PRN in the symbol.
			repeat_values<thrust::counting_iterator<int> > key_counter_avg(thrust::make_counting_iterator(0), thrust::make_counting_iterator((int)_nb_msg_prns), _prn_per_symbol);
			thrust::inclusive_scan_by_key(key_counter_avg.begin(), key_counter_avg.end(), _d_corr_mag.begin(), _d_corr_mag_avgsum.begin());

            // copy results to host
            thrust::copy(_d_corr_mag_avgsum.begin(), _d_corr_mag_avgsum.end(), _h_corr_mag_avgsum.begin());

            // search PRN index with largest magnitude and its magnitude symbol index for each PRN index in symbol
            // compute sum of all PRNs magnitudes for PRN index in symbol

            _max_avg.assign(_prn_per_symbol, 0.0);
            _mag_avgsum_sums.assign(_prn_per_symbol, 0.0);
            _max_avg_index.assign(_prn_per_symbol, 0);

            for (unsigned int i=0; i<_nb_msg_prns*_prn_per_symbol; i++)
            {
            	_mag_avgsum_sums[i%_prn_per_symbol] += _h_corr_mag_avgsum[i];

            	if (_h_corr_mag_avgsum[i] > _max_avg[i%_prn_per_symbol])
            	{
            		_max_avg[i%_prn_per_symbol] = _h_corr_mag_avgsum[i];
            		_max_avg_index[i%_prn_per_symbol] = (i/_prn_per_symbol)%_nb_msg_prns;
            	}
            }

			// loop through PRNs in symbol data to fill correlation records for this symbol

            for (unsigned int pci = 0; pci < _prn_per_symbol; pci++)
            {
                CorrelationRecord& correlation_record = pilot_correlation_analyzer.new_message_correlation_record(_transient_corr_values[pci].global_prn_i);
                
                correlation_record.selected = _transient_corr_values[pci].selected;
                correlation_record.prn_per_symbol_index = pci % _prn_per_symbol;
                correlation_record.prn_index_max = _max_avg_index[pci];
                correlation_record.magnitude_max = _max_avg[pci];
                correlation_record.noise_avg = _h_corr_mag_avgsum[(_nb_msg_prns-1)*_prn_per_symbol + pci];
                correlation_record.magnitude_avg = (_mag_avgsum_sums[pci] - correlation_record.noise_avg) / (_nb_msg_prns-1);
                correlation_record.pilot_shift = _transient_corr_values[pci].delta_t;
                correlation_record.shift_index_max = _transient_corr_values[pci].delta_t; // copy over
                correlation_record.f_rx = -_transient_corr_values[pci].delta_f; // copy over
                correlation_record.phase_at_max = _transient_corr_values[pci].pilot_phase_at_max; // copy over
            }
    	}
    }
}
