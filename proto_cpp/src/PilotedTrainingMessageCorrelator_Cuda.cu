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

     PilotedTrainingMessageCorrelator_Cuda

     Given the frequency and time displacement of correlation peak given by the
     Pilot Correlator it searches correlation for all PRNs of the training sequence symbols
     it shifts accumulates the result to detect a peak corresponding to the PRN at the 
     start of the received sequence. It uses straightforward time correlation.

     This is the CUDA implementation
     
*/

#include "PilotedTrainingMessageCorrelator_Cuda.h"
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

const cuComplex PilotedTrainingMessageCorrelator_Cuda::_c_zero = {0.0, 0.0};

PilotedTrainingMessageCorrelator_Cuda::PilotedTrainingMessageCorrelator_Cuda(
		LocalCodes_Cuda& local_codes,
		wsgc_float f_sampling,
		wsgc_float f_chip,
        unsigned int analysis_window_size,
		unsigned int sequence_length,
		unsigned int cuda_device) :
	CudaDeviceManager::CudaDeviceManager(cuda_device),
	PilotedTrainingMessageCorrelator(f_sampling, f_chip, sequence_length),
	_local_codes(local_codes),
    _local_oscillator(f_sampling, _local_codes.get_gc_generator().get_nb_code_samples(f_sampling, f_chip)),
    _nb_msg_prns(_local_codes.get_gc_generator().get_nb_message_codes()+1), // +1 for noise PRN
	_fft_N(_local_codes.get_gc_generator().get_nb_code_samples(f_sampling,f_chip)),
    _analysis_window_size(analysis_window_size),
    _max_avg(0.0),
    _max_avg_index(0),
	_d_corr_in(_fft_N),
	_d_mul_out(_fft_N*_sequence_length*analysis_window_size),
    _d_corr_out(_sequence_length*analysis_window_size),
    _d_corr_mag(_sequence_length*analysis_window_size),
    _d_corr_mag_avgsum(_sequence_length),
    _h_corr_mag_avgsum(_sequence_length),
	_d_keys(_sequence_length),
	_h_keys(_sequence_length)
{
	std::cout << "_sequence_length = " << _sequence_length << std::endl;
	assert (_sequence_length <= _nb_msg_prns);

	thrust::fill(_d_corr_out.begin(), _d_corr_out.end(), _c_zero);
	// Allocate memory areas
    _src = new wsgc_complex[_fft_N];
}


PilotedTrainingMessageCorrelator_Cuda::~PilotedTrainingMessageCorrelator_Cuda()
{
	// Free memory areas
    delete[] _src;
}


void PilotedTrainingMessageCorrelator_Cuda::execute(PilotCorrelationAnalyzer& pilot_correlation_analyzer)
{
    unsigned int pi = pilot_correlation_analyzer.get_start_analysis_pilot_correlation_records_index();
    const std::vector<PilotCorrelationRecord> pilot_correlation_records = pilot_correlation_analyzer.get_pilot_correlation_records();
    unsigned int first_training_correlation_record_index = 0;

    for (unsigned int pai=0; pai < _analysis_window_size; pi++, pai++) // scan PRNs in the window
    {
    	//CorrelationRecord& correlation_record = pilot_correlation_analyzer.new_message_correlation_record(pilot_correlation_records[pi].prn_index);
    	float delta_f = pilot_correlation_records[pi].delta_f;
    	unsigned int delta_t = pilot_correlation_records[pi].t_index_max;
        
        TrainingCorrelationRecord& correlation_record = pilot_correlation_analyzer.new_training_correlation_record(pilot_correlation_records[pi].prn_index, _sequence_length, _analysis_window_size);
        
        if (pai == 0)
        {
            first_training_correlation_record_index = pilot_correlation_analyzer.get_training_correlation_records().size() - 1;
        }
        
        correlation_record._selected = pilot_correlation_records[pi].selected;
        correlation_record._pilot_shift = delta_t;

    	// make input samples if necessary
    	if (pilot_correlation_records[pi].selected)
    	{
            std::complex<float> *prn_src = pilot_correlation_analyzer.get_synchronized_samples(pilot_correlation_records[pi].prn_index, delta_t); // get correlation peak synchronized data
        	
            if (prn_src)
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

		    	// samples member to member multiply by local codes
		    	repeat_range<thrust::device_vector<cuComplex>::iterator> d_mul_src_in(_d_corr_in.begin() ,_d_corr_in.end(), _sequence_length-pai);
		    	const thrust::device_vector<cuComplex>& d_lc_matrix = _local_codes.get_local_codes();
		    	thrust::transform(d_lc_matrix.begin()+(pai*_fft_N), d_lc_matrix.begin()+(_sequence_length*_fft_N), d_mul_src_in.begin(), _d_mul_out.begin(), cmulc_functor2());
		    	thrust::fill(_d_mul_out.begin()+((_sequence_length-pai)*_fft_N), _d_mul_out.end(), _c_zero);

		    	// sum on each PRN interval to get correlation
		    	repeat_values<thrust::counting_iterator<int> > key_sum(thrust::make_counting_iterator(0), thrust::make_counting_iterator((int)_sequence_length), _fft_N);
		    	strided_shifted_range<thrust::device_vector<cuComplex>::iterator> d_corr_out_stride(_d_corr_out.begin(), _d_corr_out.end(), _analysis_window_size, pai);
		    	thrust::reduce_by_key(key_sum.begin(), key_sum.end(), _d_mul_out.begin(), _d_keys.begin(), d_corr_out_stride.begin(), thrust::equal_to<int>(), caddc_functor());
            }
            else
            {
            	strided_shifted_range<thrust::device_vector<cuComplex>::iterator> d_corr_out_stride(_d_corr_out.begin(), _d_corr_out.end(), _analysis_window_size, pai);
            	thrust::fill(d_corr_out_stride.begin(), d_corr_out_stride.end(), _c_zero);
            }
    	}
    	else
    	{
        	strided_shifted_range<thrust::device_vector<cuComplex>::iterator> d_corr_out_stride(_d_corr_out.begin(), _d_corr_out.end(), _analysis_window_size, pai);
        	thrust::fill(d_corr_out_stride.begin(), d_corr_out_stride.end(), _c_zero);
    	}


    	//std::cout << "pai = " << pai << " pi = " << pi << std::endl;

    	// calculate averaging sums etc... for whole window at once

    	if (pai == _analysis_window_size-1)
    	{
            // calculate magnitudes (squared norm)
            thrust::transform(_d_corr_out.begin(), _d_corr_out.end(), _d_corr_mag.begin(), mag_squared_functor<cuComplex, float>());
            
			// averaging sum (that is inclusive scan with repeating iterator) of magnitudes for each PRN.
			repeat_values<thrust::counting_iterator<int> > key_counter_avg(thrust::make_counting_iterator(0), thrust::make_counting_iterator((int)_sequence_length), _analysis_window_size);
			//thrust::inclusive_scan_by_key(key_counter_avg.begin(), key_counter_avg.end(), _d_corr_mag.begin(), _d_corr_mag_avgsum.begin());
			thrust::reduce_by_key(key_counter_avg.begin(), key_counter_avg.end(), _d_corr_mag.begin(), _d_keys.begin(), _d_corr_mag_avgsum.begin());

            // copy results to host
            thrust::copy(_d_corr_mag_avgsum.begin(), _d_corr_mag_avgsum.end(), _h_corr_mag_avgsum.begin());

            // search PRN index with largest magnitude and its magnitude symbol index for each PRN index in symbol
            // compute sum of all PRNs magnitudes for PRN index in symbol

            // Prod level:
            _max_avg = 0.0;
            wsgc_float mag_avgsums_sum = 0.0;
            _max_avg_index = 0;

            for (unsigned int pi=0; pi<_sequence_length; pi++)
            {
            	//std::cout << "pi=" << pi << ": " << _h_corr_mag_avgsum[pi] << std::endl;
            	mag_avgsums_sum += _h_corr_mag_avgsum[pi];

            	if (_h_corr_mag_avgsum[pi] > _max_avg)
            	{
            		_max_avg = _h_corr_mag_avgsum[pi];
            		_max_avg_index = pi;
            	}
            }

            correlation_record._max_selected = true;

            //*/
            
            /*
            // Prototype level:
            for (unsigned int ai=0; ai<_analysis_window_size; ai++)
            {
                TrainingCorrelationRecord& correlation_record = pilot_correlation_analyzer.get_training_correlation_records()[first_training_correlation_record_index+ai];
                _max_avg = 0.0;
                wsgc_float mag_avgsums_sum = 0.0;
                _max_avg_index = 0;

                for (unsigned int pi=0; pi<_sequence_length; pi++)
                {
                    mag_avgsums_sum += _h_corr_mag_avgsum[_analysis_window_size*pi+ai];
                    
                    if (_h_corr_mag_avgsum[_analysis_window_size*pi+ai] > _max_avg)
                    {
                    	_max_avg = _h_corr_mag_avgsum[_analysis_window_size*pi+ai];
                    	_max_avg_index = pi;
                    }
                    correlation_record._max_selected = true;
                }
            }
            */
            correlation_record._magnitude_avgsum = mag_avgsums_sum-_max_avg;
            correlation_record._magnitude_max = _max_avg;
            correlation_record._prn_index_max = _max_avg_index;

            // pass results to high level correlator
            _mag_max = _max_avg;
            _maxtoavg_max = (_max_avg / correlation_record._magnitude_avgsum) * _sequence_length;
            _prni_max = _max_avg_index;
            _prnai_max = pai;

            // TODO: do something useful with the result...
            std::cout << "Max: " << _max_avg_index << " : " << _max_avg << std::endl;
    	}
    }
}
