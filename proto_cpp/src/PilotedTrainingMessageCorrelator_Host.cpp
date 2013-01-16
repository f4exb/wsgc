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

     PilotedTrainingMessageCorrelator_Host

     Given the frequency and time displacement of correlation peak given by the
     Pilot Correlator it searches correlation for all PRNs of the training sequence symbols
     it shifts accumulates the result to detect a peak corresponding to the PRN at the 
     start of the received sequence. It uses straightforward time correlation.

     This is the host implementation

*/

#include "PilotedTrainingMessageCorrelator_Host.h"
#include "GoldCodeGenerator.h"
#include "PilotCorrelationAnalyzer.h"
#include "TrainingCorrelationRecord.h"
#include "LocalCodes_Host.h"
#include "WsgcUtils.h"
#include <cmath>
#include <cstring>

PilotedTrainingMessageCorrelator_Host::PilotedTrainingMessageCorrelator_Host(
		LocalCodes_Host& local_codes,
		wsgc_float f_sampling,
		wsgc_float f_chip,
		unsigned int sequence_length) :
	PilotedTrainingMessageCorrelator(f_sampling, f_chip, sequence_length),
	_local_codes(local_codes),
    _local_oscillator(f_sampling, _local_codes.get_gc_generator().get_nb_code_samples(f_sampling, f_chip)),
    _nb_msg_prns(_local_codes.get_gc_generator().get_nb_message_codes()),
	_fft_N(_local_codes.get_gc_generator().get_nb_code_samples(f_sampling,f_chip))
{
	assert (_sequence_length <= _nb_msg_prns);

    static const wsgc_complex c_zero(0.0,0.0);

	// Allocate memory areas
    _src = new wsgc_complex[_fft_N];
	_corr_results = new wsgc_complex[_sequence_length];
    _mag_avgsums = new wsgc_float[_sequence_length];

    for (unsigned int i=0; i<_sequence_length; i++)
    {
        _mag_avgsums[i] = 0.0;
        _corr_results[i] = c_zero;
    }
}


PilotedTrainingMessageCorrelator_Host::~PilotedTrainingMessageCorrelator_Host()
{
	// Free memory areas
    delete[] _mag_avgsums;
    delete[] _corr_results;
    delete[] _src;
}


void PilotedTrainingMessageCorrelator_Host::execute(PilotCorrelationAnalyzer& pilot_correlation_analyzer)
{
    static const wsgc_complex c_zero(0.0,0.0);
    static const wsgc_complex c_one(1.0,0.0);
    wsgc_float mag_corr;
    wsgc_float mag_avgsums_sum;
    wsgc_float mag_max;
    wsgc_float corr_max;
    unsigned int prni_max = 0;
    unsigned int prni_corr_max = 0;
    const wsgc_complex *local_code;
    unsigned int ni = 0;
    
    unsigned int pi = pilot_correlation_analyzer.get_start_analysis_pilot_correlation_records_index();
    const std::vector<PilotCorrelationRecord> pilot_correlation_records = pilot_correlation_analyzer.get_pilot_correlation_records();
    _maxtoavg_max = 0.0;

    for (unsigned int i=0; i<_sequence_length; i++)
    {
        _mag_avgsums[i] = 0.0;
    }

    for (unsigned int pai=0; pai < pilot_correlation_analyzer.get_analysis_window_size_in_prns(); pi++, pai++) // scan PRNs in the window
    {
    	wsgc_float delta_f = pilot_correlation_records[pi].delta_f;
    	unsigned int delta_t = pilot_correlation_records[pi].t_index_max;
    	wsgc_complex *prn_src;
    	TrainingCorrelationRecord& correlation_record = pilot_correlation_analyzer.new_training_correlation_record(pilot_correlation_records[pi].prn_index, _sequence_length, pilot_correlation_analyzer.get_analysis_window_size_in_prns());

    	// make input samples if necessary
    	if (pilot_correlation_records[pi].selected)
    	{
    		prn_src = pilot_correlation_analyzer.get_synchronized_samples(pilot_correlation_records[pi].prn_index, delta_t); // get correlation peak synchronized data

    		// if samples are not available skip to next PRN
    		if (prn_src ==  0)
    		{
    			continue;
    		}

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
    	}
    
    	// scan all PRNs
    	mag_avgsums_sum = 0.0;
    	mag_max = 0.0;
        corr_max = 0.0;

    	for (unsigned int prni=0; prni < _sequence_length; prni++)
    	{
    		if (pilot_correlation_records[pi].selected)
    		{
    			// do correlation with PRNi
				local_code = _local_codes.get_local_code(prni);
				_corr_results[prni] = c_zero;

				for (unsigned int i=0; i < _fft_N; i++) // do the time domain correlation
				{
					//_corr_results[prni] += _src[(delta_t + i) % _fft_N] * local_code[i];
                    _corr_results[prni] += _src[i] * local_code[i]; // delta t is 0 because of synchronized samples
				}
    		}
    		else
    		{
				_corr_results[prni] = c_zero; // no contribution
    		}

    		// shifting averaging sum of magnitudes
    		WsgcUtils::magnitude_estimation(&_corr_results[prni], &mag_corr);
    		_mag_avgsums[(prni-pai+_sequence_length)%_sequence_length] += mag_corr;
    		mag_avgsums_sum += _mag_avgsums[(prni-pai+_sequence_length)%_sequence_length];
    	}
        
        // estimate point in time of message start (message epoch) at the best (max/avg) PRN available
    	for (unsigned int prni=0; prni < _sequence_length; prni++)
    	{
    		if (_mag_avgsums[prni] > mag_max)
    		{
    			mag_max = _mag_avgsums[prni];
    			prni_max = prni;
    		}
            
            WsgcUtils::magnitude_estimation(&_corr_results[prni], &mag_corr);
            if (mag_corr > corr_max)
            {
                corr_max = mag_corr;
                prni_corr_max = prni;
            }
    	}

        //std::cout << "Corr   : " << pai << " : " << corr_max << " : " << prni_corr_max << " : " << mag_max << " : " << mag_max/mag_avgsums_sum << " : " << prni_max << " : " << std::endl;

    	if (mag_max/mag_avgsums_sum > _maxtoavg_max)
    	{
    		_maxtoavg_max = mag_max/mag_avgsums_sum;
    		_mag_max = mag_max;
    		_prni_max = prni_max;
    		_prnai_max = pai;
    		correlation_record._max_selected = true;
            std::cout << "New max: " << _prnai_max << " : " << mag_max << " : " << _maxtoavg_max << " : " << _prni_max << " : " << std::endl;
    	}

		if (pilot_correlation_records[pi].selected)
		{
			correlation_record._magnitude_avgsum = mag_avgsums_sum;
			correlation_record._magnitude_max = mag_max;
			correlation_record._prn_index_max = prni_max;
		}
		else
		{
			correlation_record._magnitude_avgsum = 0.0;
			correlation_record._magnitude_max = 0.0;
			correlation_record._prn_index_max = 0;
		}

		correlation_record._pilot_shift = delta_t;
		correlation_record._selected = pilot_correlation_records[pi].selected;
    }
}

