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

     Pilot correlation analyzer

     Stores pilot correlation records and corresponding source samples for further reference. 
     Analyzes pilot correlation records to drive message correlation. Used with piloted mode.
*/

#include "PilotCorrelationAnalyzer.h"
#include "WsgcUtils.h"
#include <assert.h>
#include <cstring>
#include <iostream>
#include <algorithm>

const int PilotCorrelationAnalyzer::_time_option = CLOCK_REALTIME;
const PilotCorrelationAnalyzer::timings_t PilotCorrelationAnalyzer::tmp_time = {{0,0},{0,0}};

//=================================================================================================
PilotCorrelationAnalyzer::PilotCorrelationAnalyzer(
        unsigned int analysis_window_size, 
        unsigned int prn_per_symbol, 
        unsigned int fft_N, 
        bool two_pilots,
        unsigned int time_index_tolerance
        ) :
    _time_analyzer(fft_N, time_index_tolerance),
    _analysis_window_size(analysis_window_size),
    _prn_per_symbol(prn_per_symbol),
    _fft_N(fft_N),
    _time_index_tolerance(time_index_tolerance), //TODO: use time index tolerance
    _sample_buffer_len(_analysis_window_size*_prn_per_symbol*2),
    _global_prn_index_bot(0),
    _buffer_prn_index_bot(0),
    _buffer_prn_index_top(0),
    _buffer_prn_count(0),
    _start_pilot_correlation_records_index(0),
    _start_message_correlation_records_index(0),
    _analysis_index(0),
    _pilot_mag_display_factor(1.0),
    _message_mag_display_factor(1.0),
    _training_mag_display_factor(1.0)
{
    _samples_ext = (wsgc_complex *) WSGC_FFTW_MALLOC((_sample_buffer_len+2)*fft_N*sizeof(wsgc_fftw_complex));
    _samples = &_samples_ext[_fft_N]; // start is one PRN ahead
    std::cout << "_analysis_window_size = " << _analysis_window_size <<std::endl;
    std::cout << "_prn_per_symbol = " << _prn_per_symbol << std::endl;
    std::cout << "_sample_buffer_len = " << _sample_buffer_len << std::endl;
}


//=================================================================================================
PilotCorrelationAnalyzer::~PilotCorrelationAnalyzer()
{
    WSGC_FFTW_FREE(_samples_ext);
}


//=================================================================================================
void PilotCorrelationAnalyzer::store_prn_samples(wsgc_complex *samples)
{
    if (_buffer_prn_count < _sample_buffer_len)
    {
        _buffer_prn_index_top = (_buffer_prn_index_bot + _buffer_prn_count) % _sample_buffer_len;
        _buffer_prn_count++;
    }
    else
    {
        _buffer_prn_index_top = _buffer_prn_index_bot;
        _buffer_prn_index_bot = (_buffer_prn_index_bot + 1) % _sample_buffer_len;
        _global_prn_index_bot++;
    }

	//std::cout << "store_prn_samples: _global_prn_index_bot=" << _global_prn_index_bot << ", _buffer_prn_index_top=" << _buffer_prn_index_top << ", _buffer_prn_count=" << _buffer_prn_count << std::endl;
	memcpy((void *) &_samples[_buffer_prn_index_top*_fft_N], (void *) samples, _fft_N*sizeof(wsgc_complex));

	if (_buffer_prn_index_top == 0) // first element
	{
		memcpy((void *) &_samples[_sample_buffer_len*_fft_N], (void *) samples, _fft_N*sizeof(wsgc_complex)); // re-copy at extra end
	}
}


//=================================================================================================
wsgc_complex *PilotCorrelationAnalyzer::get_samples(unsigned int global_prn_index)
{
    if (global_prn_index < _global_prn_index_bot) // no more in buffer
    {
    	std::cout << "get_samples: global_prn_index=" << global_prn_index << ", _global_prn_index_bot=" << _global_prn_index_bot << ", _buffer_prn_index_bot=" << _buffer_prn_index_bot << ", _buffer_prn_count=" << _buffer_prn_count << std::endl;
        std::cout << "    + not in buffer anymore" << std::endl;
        return 0;
    }
    else if (global_prn_index >= _global_prn_index_bot + _buffer_prn_count) // not yet in buffer
    {
    	std::cout << "get_samples: global_prn_index=" << global_prn_index << ", _global_prn_index_bot=" << _global_prn_index_bot << ", _buffer_prn_index_bot=" << _buffer_prn_index_bot << ", _buffer_prn_count=" << _buffer_prn_count << std::endl;
        std::cout << "    + not yet in buffer" << std::endl;
        return 0;
    }
    else
    {
        unsigned int source_index = ((global_prn_index - _global_prn_index_bot) + _buffer_prn_index_bot) % _sample_buffer_len;
    
        if (source_index == _sample_buffer_len - 1) // last element
        {
        	memcpy((void *) _samples_ext, (void *) &_samples[source_index*_fft_N], _fft_N*sizeof(wsgc_complex)); // re-copy at extra begin
        }

        //std::cout << "    + return: " << source_index << std::endl;
        
        return &_samples[source_index*_fft_N];
    }
}


//=================================================================================================
wsgc_complex *PilotCorrelationAnalyzer::get_synchronized_samples(unsigned int global_prn_index, unsigned int shift_index)
{
	bool use_preceding_samples = shift_index > _fft_N / 2;
	unsigned int source_index = (global_prn_index - _global_prn_index_bot + _buffer_prn_index_bot + (use_preceding_samples ? -1 : 0)) % _sample_buffer_len;

    if (global_prn_index < _global_prn_index_bot) // no more in buffer
    {
    	std::cout << "get_samples: global_prn_index=" << global_prn_index << ", _global_prn_index_bot=" << _global_prn_index_bot << ", _buffer_prn_index_bot=" << _buffer_prn_index_bot << ", _buffer_prn_count=" << _buffer_prn_count << ", result=" << source_index << std::endl;
        std::cout << "    + not in buffer anymore" << std::endl;
        return 0;
    }
    else if (global_prn_index >= _global_prn_index_bot + _buffer_prn_count) // not yet in buffer
    {
    	std::cout << "get_samples: global_prn_index=" << global_prn_index << ", _global_prn_index_bot=" << _global_prn_index_bot << ", _buffer_prn_index_bot=" << _buffer_prn_index_bot << ", _buffer_prn_count=" << _buffer_prn_count << ", result=" << source_index << std::endl;
        std::cout << "    + not yet in buffer" << std::endl;
        return 0;
    }
    else // current is in buffer
    {
        if ((global_prn_index - _global_prn_index_bot + _buffer_prn_index_bot) == _sample_buffer_len - 1) // last element
        {
        	memcpy((void *) _samples_ext, (void *) &_samples[source_index*_fft_N], _fft_N*sizeof(wsgc_complex)); // re-copy at extra begin
        }

        if (use_preceding_samples)
        {
            if (global_prn_index == _global_prn_index_bot)
            {
            	std::cout << "get_samples: global_prn_index=" << global_prn_index << ", _global_prn_index_bot=" << _global_prn_index_bot << ", _buffer_prn_index_bot=" << _buffer_prn_index_bot << ", _buffer_prn_count=" << _buffer_prn_count << ", result=" << source_index << std::endl;
                std::cout << "    + previous not in buffer anymore" << std::endl;
                return 0; // cannot serve before start. previous no more in buffer
            }
            else 
            {
                //std::cout << "    + return: " <<source_index << "(" << shift_index << ")" << std::endl;
                return &_samples[source_index*_fft_N + shift_index];
            }
        }
        else 
        {
            if (global_prn_index == _global_prn_index_bot + _buffer_prn_count - 1)
            {
            	std::cout << "get_samples: global_prn_index=" << global_prn_index << ", _global_prn_index_bot=" << _global_prn_index_bot << ", _buffer_prn_index_bot=" << _buffer_prn_index_bot << ", _buffer_prn_count=" << _buffer_prn_count << ", result=" << source_index << std::endl;
                std::cout << "    + next not yet in buffer" << std::endl;
                return 0; // cannot serve past end. next is not yet in buffer
            }
            else 
            {
                //std::cout << "    + return: " <<source_index << "(" << shift_index << ")" << std::endl;
                return &_samples[source_index*_fft_N + shift_index];
            }
        }
    }
}


//=================================================================================================
const wsgc_complex *PilotCorrelationAnalyzer::get_last_samples() const
{
	return &_samples[(_buffer_prn_index_top)*_fft_N];
}


//=================================================================================================
PilotCorrelationRecord& PilotCorrelationAnalyzer::new_pilot_correlation_record()
{
    return new_pilot_correlation_record(_pilot_correlation_records);
}


//=================================================================================================
CorrelationRecord& PilotCorrelationAnalyzer::new_message_correlation_record(unsigned int global_prn_index)
{
    static const CorrelationRecord tmp_message_correlation_record;

    _message_correlation_records.push_back(tmp_message_correlation_record);
    _message_correlation_records.back().global_prn_index = global_prn_index;
    _message_correlation_records.back().prn_per_symbol = _prn_per_symbol;

    return _message_correlation_records.back();
}


//=================================================================================================
TrainingCorrelationRecord& PilotCorrelationAnalyzer::new_training_correlation_record(
		unsigned int global_prn_index,
		unsigned int sequence_length,
		unsigned int analysis_window_prn)

{
    static const TrainingCorrelationRecord tmp_training_correlation_record(sequence_length, analysis_window_prn);

    _training_correlation_records.push_back(tmp_training_correlation_record);
    _training_correlation_records.back()._global_prn_index = global_prn_index;

    return _training_correlation_records.back();
}


//=================================================================================================
PilotCorrelationRecord& PilotCorrelationAnalyzer::new_pilot_correlation_record(std::vector<PilotCorrelationRecord>& pilot_correlation_records)
{
    static const PilotCorrelationRecord tmp_pilot_correlation_record;

    pilot_correlation_records.push_back(tmp_pilot_correlation_record);
    pilot_correlation_records.back().prn_index = get_prn_index();

    return pilot_correlation_records.back();
}


//=================================================================================================
PilotCorrelationRecord& PilotCorrelationAnalyzer::get_pilot_correlation_record_back(unsigned int reverse_index)
{
    return get_pilot_correlation_record_back(_pilot_correlation_records, reverse_index);
}


//=================================================================================================
PilotCorrelationRecord& PilotCorrelationAnalyzer::get_pilot_correlation_record_back(std::vector<PilotCorrelationRecord>& pilot_correlation_records, unsigned int reverse_index)
{
    return pilot_correlation_records.at(pilot_correlation_records.size()-1-reverse_index);
}


//=================================================================================================
bool PilotCorrelationAnalyzer::validate_pilot_correlation_records_back(unsigned int reverse_index)
{
    validate_pilot_correlation_record(get_pilot_correlation_record_back(_pilot_correlation_records, reverse_index));
    return _time_analyzer.get_validate_count() == _analysis_window_size * _prn_per_symbol;
}


//=================================================================================================
void PilotCorrelationAnalyzer::validate_pilot_correlation_record(
    PilotCorrelationRecord& pilot_correlation_record)
{
	_time_analyzer.validate_time_correlation_record(pilot_correlation_record);
}


//=================================================================================================
bool PilotCorrelationAnalyzer::analyze(unsigned int& best_time_shift_start, unsigned int& best_time_shift_length)
{
	return _time_analyzer.analyze_time_shifts(_pilot_correlation_records, _start_pilot_correlation_records_index, best_time_shift_start, best_time_shift_length);
}


//=================================================================================================
void PilotCorrelationAnalyzer::post_process_noise()
{
	wsgc_float noise_sum = 0.0;
	unsigned int valid_count = 0;

	for (unsigned int i = _start_message_correlation_records_index; i < _message_correlation_records.size(); i++)
	{
		if (_message_correlation_records[i].selected)
		{
			noise_sum += _message_correlation_records[i].noise_avg;
			valid_count++;
		}
	}

	for (unsigned int i = _start_message_correlation_records_index; i < _message_correlation_records.size(); i++)
	{
		_message_correlation_records[i].noise_avg = noise_sum / valid_count;
	}
}


//=================================================================================================
void PilotCorrelationAnalyzer::reset_analysis()
{
    _time_analyzer.reset_analysis();
    _start_pilot_correlation_records_index = (_pilot_correlation_records.size() / (_analysis_window_size*_prn_per_symbol)) * _analysis_window_size * _prn_per_symbol; // record beginning of new analysis window in the pilot correlation records
    _start_message_correlation_records_index = _message_correlation_records.size(); // record beginning of new analysis window in the message correlation records
    _analysis_index++;
}


//=================================================================================================
void PilotCorrelationAnalyzer::dump_pilot_correlation_records(std::ostringstream& os) const
{
	std::vector<PilotCorrelationRecord>::const_iterator it;
	std::vector<PilotCorrelationRecord>::const_iterator it_end;

	it = _pilot_correlation_records.begin();
	it_end = _pilot_correlation_records.end();

	PilotCorrelationRecord::dump_oneline_banner(os);

	for (; it != it_end; ++it)
	{
		it->dump_oneline(os, _pilot_mag_display_factor);
	}
}


//=================================================================================================
void PilotCorrelationAnalyzer::dump_message_correlation_records(std::ostringstream& os) const
{
	std::vector<CorrelationRecord>::const_iterator it = _message_correlation_records.begin();
	std::vector<CorrelationRecord>::const_iterator it_end = _message_correlation_records.end();

	CorrelationRecord::dump_banner(os);

	for (; it != it_end; ++it)
	{
		it->dump_line(_message_mag_display_factor, os);
	}
}


//=================================================================================================
void PilotCorrelationAnalyzer::dump_training_correlation_records(std::ostringstream& os) const
{
	std::vector<TrainingCorrelationRecord>::const_iterator it = _training_correlation_records.begin();
	std::vector<TrainingCorrelationRecord>::const_iterator it_end = _training_correlation_records.end();

	TrainingCorrelationRecord::dump_banner(os);

	for (; it != it_end; ++it)
	{
		it->dump_line(os, _training_mag_display_factor);
	}
}


//=================================================================================================
void PilotCorrelationAnalyzer::dump_histo_time_shift_occurences(std::ostringstream& os) const
{
	_time_analyzer.dump_histo_time_shift_occurences(os);
}


//=================================================================================================
void PilotCorrelationAnalyzer::dump_timings(std::ostringstream& os)
{
	//----------------------//
	os << "Pilot multiplication and IFFT timings:" << std::endl;
	os << "First->Last : " << std::setw(12) << std::setprecision(9) << WsgcUtils::get_time_difference(_pilot_mul_ifft_times.back()._end, _pilot_mul_ifft_times.front()._beg) << " s" << std::endl;
	std::vector<timings_t>::const_iterator pit = _pilot_mul_ifft_times.begin();
	const std::vector<timings_t>::const_iterator pend = _pilot_mul_ifft_times.end();
	unsigned int i=0;
	double timing_sum = 0.0;
	double timing;

	for (; pit < pend; ++pit, i++)
	{
		timing = WsgcUtils::get_time_difference(pit->_end, pit->_beg);
		timing_sum += timing;
		os << std::setw(3) <<  i << " : " << std::left << std::setw(12) << std::setprecision(9) << timing << " s" << std::endl;
	}

	os << " Total : " << std::setw(12) << std::setprecision(9) << timing_sum << " s" << std::endl;

	//----------------------//
	os << "Pilot averaging and reduction timings:" << std::endl;
	os << "First->Last : " << std::setw(12) << std::setprecision(9) << WsgcUtils::get_time_difference(_pilot_avg_times.back()._end, _pilot_avg_times.front()._beg) << " s" << std::endl;
	std::vector<timings_t>::const_iterator ait = _pilot_avg_times.begin();
	const std::vector<timings_t>::const_iterator aend = _pilot_avg_times.end();
	i=0;
	timing_sum = 0.0;

	for (; ait < aend; ++ait, i++)
	{
		timing = WsgcUtils::get_time_difference(ait->_end, ait->_beg);
		timing_sum += timing;
		os << std::setw(3) <<  i << " : " << std::left << std::setw(12) << std::setprecision(9) << timing << " s" << std::endl;
	}

	os << " Total : " << std::setw(12) << std::setprecision(9) << timing_sum << " s" << std::endl;

	os << "Message timings:" << std::endl;

	//----------------------//
	os << "Pilot CUDA averaging timings:" << std::endl;
	std::vector<timings_t>::const_iterator cait = _pilot_cuda_avg_times.begin();
	const std::vector<timings_t>::const_iterator caend = _pilot_cuda_avg_times.end();
	i=0;
	timing_sum = 0.0;

	for (; cait < caend; ++cait, i++)
	{
		timing = WsgcUtils::get_time_difference(cait->_end, cait->_beg);
		timing_sum += timing;
		os << std::setw(3) <<  i << " : " << std::left << std::setw(12) << std::setprecision(9) << timing << " s" << std::endl;
	}

	os << " Total : " << std::setw(12) << std::setprecision(9) << timing_sum << " s" << std::endl;

	//----------------------//
	os << "Pilot CUDA reduction timings:" << std::endl;
	std::vector<timings_t>::const_iterator crit = _pilot_cuda_reduc_times.begin();
	const std::vector<timings_t>::const_iterator crend = _pilot_cuda_reduc_times.end();
	i=0;
	timing_sum = 0.0;

	for (; crit < crend; ++crit, i++)
	{
		timing = WsgcUtils::get_time_difference(crit->_end, crit->_beg);
		timing_sum += timing;
		os << std::setw(3) <<  i << " : " << std::left << std::setw(12) << std::setprecision(9) << timing << " s" << std::endl;
	}

	os << " Total : " << std::setw(12) << std::setprecision(9) << timing_sum << " s" << std::endl;

	//----------------------//
	os << "Message timings:" << std::endl;

	std::vector<timings_t>::const_iterator mit = _message_times.begin();
	const std::vector<timings_t>::const_iterator mend = _message_times.end();
	i=0;

	for (; mit < mend; ++mit, i++)
	{
		os << std::setw(3) <<  i << " : " << std::left << std::setw(12) << std::setprecision(9) << WsgcUtils::get_time_difference(mit->_end, mit->_beg) << " s" << std::endl;
	}
}
