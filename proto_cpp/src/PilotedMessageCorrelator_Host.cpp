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

#include "PilotedMessageCorrelator_Host.h"
#include "GoldCodeGenerator.h"
#include "CorrelationRecord.h"
#include "PilotCorrelationAnalyzer.h"
#include "LocalCodes_Host.h"
#include "WsgcUtils.h"
#include <cmath>
#include <cstring>

PilotedMessageCorrelator_Host::PilotedMessageCorrelator_Host(
		LocalCodes_Host& local_codes,
		wsgc_float f_sampling,
		wsgc_float f_chip,
		unsigned int prn_per_symbol) :
	PilotedMessageCorrelator(f_sampling, f_chip, prn_per_symbol),
	_local_codes(local_codes),
    _local_oscillator(f_sampling, _local_codes.get_gc_generator().get_nb_code_samples(f_sampling, f_chip)),
    _nb_msg_prns(_local_codes.get_gc_generator().get_nb_message_codes()+1), // +1 for noise PRN
	_fft_N(_local_codes.get_gc_generator().get_nb_code_samples(f_sampling,f_chip))
{
    static const wsgc_complex c_zero(0.0,0.0);

	// Allocate memory areas
    _src = new wsgc_complex[_fft_N];
	_corr_results = new wsgc_complex[_nb_msg_prns*_prn_per_symbol];
    _corr_avgsums = new wsgc_complex[_nb_msg_prns];
	_noise_corr_results = new wsgc_complex[_prn_per_symbol];

    for (unsigned int i=0; i<_nb_msg_prns; i++)
    {
        _corr_avgsums[i] = c_zero;
    
        for (unsigned int j=0; j<_prn_per_symbol; j++)
        {
            _corr_results[i*_prn_per_symbol + j] = c_zero;
        }
    }

    for (unsigned int j=0; j<_prn_per_symbol; j++)
    {
    	_noise_corr_results[j] = c_zero;
    }
}


PilotedMessageCorrelator_Host::~PilotedMessageCorrelator_Host()
{
	// Free memory areas
	delete[] _noise_corr_results;
    delete[] _corr_avgsums;
    delete[] _corr_results;
    delete[] _src;
}


void PilotedMessageCorrelator_Host::execute(PilotCorrelationAnalyzer& pilot_correlation_analyzer)
{
    static const wsgc_complex c_zero(0.0,0.0);
    static const wsgc_complex c_one(1.0,0.0);
    wsgc_complex corr_result;
    wsgc_float max_magnitude, max_magnitude_i;
    wsgc_float avgsum_magnitude;
    wsgc_complex max_correlation;
    wsgc_complex noise_avgsum;
    wsgc_float magnitude, magnitude_sum, noise_magnitude, avg_noise;
    unsigned int max_prn_index, max_prn_index_i;
    const wsgc_complex *local_code;
    unsigned int ni = 0;
    
    unsigned int pi = pilot_correlation_analyzer.get_start_analysis_pilot_correlation_records_index();
    const std::vector<PilotCorrelationRecord> pilot_correlation_records = pilot_correlation_analyzer.get_pilot_correlation_records();

    for (unsigned int pai=0; pai < pilot_correlation_analyzer.get_analysis_window_size()*_prn_per_symbol; pi++, pai++) // scan PRNs in the window
    {
    	unsigned int ai = pai % _prn_per_symbol;
    	CorrelationRecord& correlation_record = pilot_correlation_analyzer.new_message_correlation_record(pilot_correlation_records[pi].prn_index);
    	wsgc_float delta_f = pilot_correlation_records[pi].delta_f;
    	unsigned int delta_t = pilot_correlation_records[pi].t_index_max;
    	wsgc_complex *prn_src = pilot_correlation_analyzer.get_samples(pilot_correlation_records[pi].prn_index);

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

    	max_magnitude = 0.0;
    	max_magnitude_i = 0.0;
    	avgsum_magnitude = 0.0;
    	max_prn_index = _nb_msg_prns-1;
    	max_correlation = c_one;

    	// Zero sum averaging memory on averaging start
    	if (ai==0)
    	{
    		for (unsigned int prni=0; prni < _nb_msg_prns; prni++)
    		{
    			for (unsigned int aj=0; aj<_prn_per_symbol; aj++)
    			{
    				_corr_results[prni*_prn_per_symbol + aj] = c_zero;
    			}
    		}
    	}

    	// scan all PRNs
    	for (unsigned int prni=0; prni < _nb_msg_prns; prni++)
    	{
    		if (pilot_correlation_records[pi].selected)
    		{
    			// do correlation with PRNi
				local_code = _local_codes.get_local_code(prni);
				_corr_results[prni*_prn_per_symbol + ai] = c_zero;

				for (unsigned int i=0; i < _fft_N; i++) // do the time domain correlation
				{
					_corr_results[prni*_prn_per_symbol + ai] += _src[(delta_t + i) % _fft_N] * local_code[i];
				}

				if (prni == _nb_msg_prns-1) // Noise PRN - copy result in noise area, only non null results contribute
				{
					memcpy(&_noise_corr_results[ni], &_corr_results[prni*_prn_per_symbol + ai], sizeof(wsgc_complex));
					ni = (ni+1) % _prn_per_symbol;
				}
    		}
    		else
    		{
				_corr_results[prni*_prn_per_symbol + ai] = c_zero; // no contribution
    		}

    		// instant and averaging calculations
			if (prni < _nb_msg_prns-1) // real message PRNs
			{
				/*// vector sum
				_corr_avgsums[prni] = _corr_results[prni*_prn_per_symbol];

				for (unsigned int aj=1; aj<_prn_per_symbol; aj++)
				{
					_corr_avgsums[prni] += _corr_results[prni*_prn_per_symbol + aj];
				}

				// maximum magnitude and average
				WsgcUtils::magnitude_estimation(&_corr_avgsums[prni], &magnitude);
				*/

				// instant magnitude
				WsgcUtils::magnitude_estimation(&_corr_results[prni*_prn_per_symbol + ai], &magnitude);

				if (magnitude == 0)
				{
					max_prn_index_i = _nb_msg_prns-1; // not selected => set as noise PRN
				}
				else if (magnitude > max_magnitude_i)
				{
					max_magnitude_i = magnitude;
					max_prn_index_i = prni;
				}

				// magnitude sum
				magnitude_sum = 0.0;
				for (unsigned int aj=0; aj<_prn_per_symbol; aj++)
				{
					WsgcUtils::magnitude_estimation(&_corr_results[prni*_prn_per_symbol + aj], &magnitude);
					magnitude_sum += magnitude;
				}
				magnitude = magnitude_sum;

				avgsum_magnitude += magnitude;

				if (magnitude > max_magnitude)
				{
					max_magnitude = magnitude;
					max_prn_index = prni;
					max_correlation = _corr_avgsums[prni];
				}
			}
			else // noise PRN
			{
				avg_noise = 0.0;

				for (unsigned int aj=0; aj<_prn_per_symbol; aj++)
				{
					WsgcUtils::magnitude_estimation(&_noise_corr_results[aj], &noise_magnitude);
					avg_noise += noise_magnitude;
				}
			}
    	}

    	correlation_record.prn_per_symbol_index = ai;
    	correlation_record.prn_index_max_i = max_prn_index_i;
    	correlation_record.magnitude_max_i = max_magnitude_i;
    	correlation_record.prn_index_max = max_prn_index;
    	correlation_record.magnitude_max = max_magnitude;
    	correlation_record.magnitude_avg = avgsum_magnitude / (_nb_msg_prns-1);
    	correlation_record.noise_avg = avg_noise;
    	correlation_record.pilot_shift = delta_t;
    	correlation_record.shift_index_max = delta_t; // copy over
    	correlation_record.f_rx = -delta_f; // copy over
    	correlation_record.phase_at_max = atan2(max_correlation.imag(), max_correlation.real());
    }
}

