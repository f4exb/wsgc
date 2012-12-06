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
const wsgc_float PilotCorrelationAnalyzer::_best_time_margin_threshold = 0.13; // selected correlation time delta difference with next / number of correlations ratio threshold

//=================================================================================================
PilotCorrelationAnalyzer::PilotCorrelationAnalyzer(
        unsigned int analysis_window_size, 
        unsigned int prn_per_symbol, 
        unsigned int fft_N, 
        bool two_pilots,
        unsigned int time_index_tolerance
        ) :
    _analysis_window_size(analysis_window_size),
    _prn_per_symbol(prn_per_symbol),
    _fft_N(fft_N),
    _two_pilots(two_pilots),
    _time_index_tolerance(time_index_tolerance),
    _sample_buffer_len(2*analysis_window_size*prn_per_symbol),
    _samples_start_global_prn_index(0),
    _prn_index(0),
    _start_pilot_correlation_records_index(0),
    _sum_max_pilot1(0),
    _sum_max_pilot2(0),
    _start_message_correlation_records_index(0),
    _validate_count(0),
    _analysis_index(0)
{
    _samples = (wsgc_complex *) WSGC_FFTW_MALLOC(analysis_window_size*prn_per_symbol*fft_N*2*sizeof(wsgc_fftw_complex));
}


//=================================================================================================
PilotCorrelationAnalyzer::~PilotCorrelationAnalyzer()
{
    WSGC_FFTW_FREE(_samples);
}


//=================================================================================================
void PilotCorrelationAnalyzer::store_prn_samples(wsgc_complex *samples)
{
	if (_prn_index >= _sample_buffer_len)
	{
		_prn_index = 0;
		_samples_start_global_prn_index += _sample_buffer_len;

	}

	memcpy((void *) &_samples[_prn_index*_fft_N], (void *) samples, _fft_N*sizeof(wsgc_complex));
	_prn_index++;
}


//=================================================================================================
wsgc_complex *PilotCorrelationAnalyzer::get_samples(unsigned int global_prn_index)
{
	if (global_prn_index < _samples_start_global_prn_index)
	{
		if (global_prn_index < _samples_start_global_prn_index - _sample_buffer_len + _prn_index)
		{
			return 0; // gone
		}
		else
		{
			return &_samples[(global_prn_index - _samples_start_global_prn_index + _sample_buffer_len)*_fft_N];
		}
	}
	else
	{
		if ((global_prn_index - _samples_start_global_prn_index) < _prn_index)
		{
			return &_samples[(global_prn_index - _samples_start_global_prn_index)*_fft_N];
		}
		else
		{
			return 0; // not there yet
		}
	}

	return 0;
}


//=================================================================================================
wsgc_complex *PilotCorrelationAnalyzer::get_last_samples()
{
	return &_samples[(_prn_index - 1)*_fft_N];
}


//=================================================================================================
PilotCorrelationRecord& PilotCorrelationAnalyzer::new_pilot_correlation_record(bool alternate)
{
    if (alternate)
    {
        return new_pilot_correlation_record(_pilot2_correlation_records);
    }
    else
    {
        return new_pilot_correlation_record(_pilot1_correlation_records);
    }
}


//=================================================================================================
CorrelationRecord& PilotCorrelationAnalyzer::new_message_correlation_record(unsigned int global_prn_index)
{
    static const CorrelationRecord tmp_message_correlation_record;

    _message_correlation_records.push_back(tmp_message_correlation_record);
    _message_correlation_records.back().global_prn_index = global_prn_index;

    return _message_correlation_records.back();
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
PilotCorrelationRecord& PilotCorrelationAnalyzer::get_pilot_correlation_record_back(unsigned int reverse_index, bool alternate)
{
    if (alternate)
    {
        return get_pilot_correlation_record_back(_pilot2_correlation_records, reverse_index);
    }
    else
    {
        return get_pilot_correlation_record_back(_pilot1_correlation_records, reverse_index);
    }
}


//=================================================================================================
PilotCorrelationRecord& PilotCorrelationAnalyzer::get_pilot_correlation_record_back(std::vector<PilotCorrelationRecord>& pilot_correlation_records, unsigned int reverse_index)
{
    return pilot_correlation_records.at(pilot_correlation_records.size()-1-reverse_index);
}


//=================================================================================================
bool PilotCorrelationAnalyzer::validate_pilot_correlation_records_back(unsigned int reverse_index)
{
    assert((!_two_pilots) || (_pilot1_correlation_records.size() == _pilot2_correlation_records.size())); // ensure both pilot records are at the same point if both are present

    validate_pilot_correlation_record(get_pilot_correlation_record_back(_pilot1_correlation_records, reverse_index), _pilot1_time_shift_occurences, _sum_max_pilot1);
    _validate_count++;

    if (_two_pilots)
    {
        validate_pilot_correlation_record(get_pilot_correlation_record_back(_pilot2_correlation_records, reverse_index), _pilot2_time_shift_occurences, _sum_max_pilot2);
    }
    
    return _validate_count == _analysis_window_size * _prn_per_symbol;
}


//=================================================================================================
void PilotCorrelationAnalyzer::validate_pilot_correlation_record(
    PilotCorrelationRecord& pilot_correlation_record, 
    std::map<unsigned int, unsigned int>& shift_occurences_dict,
    wsgc_float& sum_max_pilot)
{
    // store current shift maximum in shifts dictionnary
    unsigned int shift_index_max = pilot_correlation_record.t_index_max;
    std::map<unsigned int, unsigned int>::iterator shift_occurences_it;
    shift_occurences_it = shift_occurences_dict.find(shift_index_max);
    
    if (shift_occurences_it == shift_occurences_dict.end())
    {
        shift_occurences_dict[shift_index_max] = 1;
    }
    else
    {
        shift_occurences_dict[shift_index_max] += 1;
    }

    // accumulate the sum of maximum magnitudes
    sum_max_pilot += pilot_correlation_record.magnitude_max;
}


//=================================================================================================
bool PilotCorrelationAnalyzer::analyze(bool& alternate_pilot, unsigned int& best_time_shift_start, unsigned int& best_time_shift_length)
{
    assert((!_two_pilots) || (_pilot1_correlation_records.size() == _pilot2_correlation_records.size())); // ensure both pilot records are at the same point if both are present

    if ((!_two_pilots) || (_sum_max_pilot1 > _sum_max_pilot2)) // pilot1 selected
    {
        alternate_pilot = false;
        _pilot_numbers.push_back(1);
        return analyze_time_shifts(_pilot1_time_shift_occurences, _pilot1_correlation_records, best_time_shift_start, best_time_shift_length);
    }
    else
    {
        alternate_pilot = true;
        _pilot_numbers.push_back(2);
        return analyze_time_shifts(_pilot2_time_shift_occurences, _pilot2_correlation_records, best_time_shift_start, best_time_shift_length);
    }
}


//=================================================================================================
void PilotCorrelationAnalyzer::make_time_shift_histogram(std::map<unsigned int, unsigned int>& shift_occurences_dict)
{
	static const std::vector<PrnTimeShiftRecord> tmp_histo_time_shift_occurences;
    _histo_time_shift_occurences_vector.push_back(tmp_histo_time_shift_occurences);

    std::map<unsigned int, unsigned int>::const_iterator shift_occurences_it = shift_occurences_dict.begin();
    const std::map<unsigned int, unsigned int>::const_iterator shift_occurences_end = shift_occurences_dict.end();
    std::map<unsigned int, unsigned int>::const_iterator last_shift_occurence_it = shift_occurences_dict.begin();
    static const PrnTimeShiftRecord tmp_time_shift_record;
    static const std::vector<PrnTimeShiftRecord>::iterator _histo_end = _histo_time_shift_occurences_vector.back().end();
    bool zero_shift = false;
    unsigned int zero_shift_index, last_shift_index;
    bool last_shift = false;

    for (; shift_occurences_it != shift_occurences_end; ++shift_occurences_it)
    {
        if (shift_occurences_it->first == ((last_shift_occurence_it->first) + 1) % _fft_N)
        {
            _histo_time_shift_occurences_vector.back().back().peak_mag_sum += shift_occurences_it->second;
            _histo_time_shift_occurences_vector.back().back().peak_length++;
            
            if (shift_occurences_it->second > _histo_time_shift_occurences_vector.back().back().peak_max)
            {
            	_histo_time_shift_occurences_vector.back().back().peak_max = shift_occurences_it->second;
            	_histo_time_shift_occurences_vector.back().back().peak_shift_at_max = shift_occurences_it->first;
            }
        }
        else
        {
        	_histo_time_shift_occurences_vector.back().push_back(tmp_time_shift_record);
        	_histo_time_shift_occurences_vector.back().back().peak_start = shift_occurences_it->first;
        	_histo_time_shift_occurences_vector.back().back().peak_mag_sum = shift_occurences_it->second;
        	_histo_time_shift_occurences_vector.back().back().peak_length = 1;
        	_histo_time_shift_occurences_vector.back().back().peak_shift_at_max = shift_occurences_it->first;
        	_histo_time_shift_occurences_vector.back().back().peak_max = shift_occurences_it->second;
        }
        
        if (shift_occurences_it->first == 0) // potentially combines with item at end of sequence
        {
            zero_shift = true;
            zero_shift_index = _histo_time_shift_occurences_vector.back().size() - 1;
        }

        if (shift_occurences_it->first == _fft_N - 1) // potentially combines with item at start of sequence
        {
            last_shift = true;
            last_shift_index = _histo_time_shift_occurences_vector.back().size() - 1;
        }
            
        last_shift_occurence_it = shift_occurences_it;
    }
    
    // coalesce intervals at extremities of the time shift window
    if (zero_shift && last_shift)
    {
    	_histo_time_shift_occurences_vector.back()[last_shift_index].peak_mag_sum += _histo_time_shift_occurences_vector.back()[zero_shift_index].peak_mag_sum;
    	_histo_time_shift_occurences_vector.back()[last_shift_index].peak_length += _histo_time_shift_occurences_vector.back()[zero_shift_index].peak_length;

    	if (_histo_time_shift_occurences_vector.back()[zero_shift_index].peak_max > _histo_time_shift_occurences_vector.back()[last_shift_index].peak_max)
    	{
    		_histo_time_shift_occurences_vector.back()[last_shift_index].peak_max = _histo_time_shift_occurences_vector.back()[zero_shift_index].peak_max;
    		_histo_time_shift_occurences_vector.back()[last_shift_index].peak_shift_at_max = _histo_time_shift_occurences_vector.back()[zero_shift_index].peak_shift_at_max;
    	}

    	_histo_time_shift_occurences_vector.back().erase(_histo_time_shift_occurences_vector.back().begin() + zero_shift_index);
    }
    
    std::sort(_histo_time_shift_occurences_vector.back().begin(), _histo_time_shift_occurences_vector.back().end(), PrnTimeShiftRecord::order_peak_mag_sum);
}


//=================================================================================================
bool PilotCorrelationAnalyzer::analyze_time_shifts(
    std::map<unsigned int, unsigned int>& shift_occurences_dict,
    std::vector<PilotCorrelationRecord>& pilot_correlation_records, 
    unsigned int& best_time_shift_start,
    unsigned int& best_time_shift_length)
{
    make_time_shift_histogram(shift_occurences_dict);
    
    PrnTimeShiftRecord& best_time_shift_peak = _histo_time_shift_occurences_vector.back().front();
    best_time_shift_start = best_time_shift_peak.peak_start;
    best_time_shift_length = best_time_shift_peak.peak_length;
    unsigned int best_time_shift_margin;

    if (_histo_time_shift_occurences_vector.back().size() > 1)
    {
    	best_time_shift_margin = best_time_shift_peak.peak_mag_sum - (_histo_time_shift_occurences_vector.back().begin()+1)->peak_mag_sum;
    }
    else
    {
    	best_time_shift_margin = best_time_shift_peak.peak_mag_sum;
    }

    bool time_shift_margin_ko = ((wsgc_float) best_time_shift_margin/_validate_count) < _best_time_margin_threshold;

    std::vector<PilotCorrelationRecord>::iterator correlation_record_it = pilot_correlation_records.begin() + _start_pilot_correlation_records_index;
    const std::vector<PilotCorrelationRecord>::iterator correlation_record_end = pilot_correlation_records.end();
      
    for (; correlation_record_it != correlation_record_end; ++correlation_record_it)
    {
    	if (time_shift_margin_ko)
    	{
    		correlation_record_it->selected = false;
    	}
    	else
    	{
            // checks whether the PRN time shift lies in the best PRN time shift peak interval
            if (best_time_shift_peak.peak_start + best_time_shift_peak.peak_length < _fft_N) // no time shift window boundary in the middle
            {
                if ((correlation_record_it->t_index_max < best_time_shift_peak.peak_start) 
                 || (correlation_record_it->t_index_max >= best_time_shift_peak.peak_start + best_time_shift_peak.peak_length))
                {
                    correlation_record_it->selected = false;
                }
                else
                {
                    correlation_record_it->selected = true;
                }
            }
            else
            {
                if ((correlation_record_it->t_index_max < best_time_shift_peak.peak_start) 
                 && (correlation_record_it->t_index_max >= (best_time_shift_peak.peak_start + best_time_shift_peak.peak_length) % _fft_N))
                {
                    correlation_record_it->selected = false;
                }
                else
                {
                    correlation_record_it->selected = true;
                }
            }
    	}
    }
    
    return !time_shift_margin_ko;
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
    assert((!_two_pilots) || (_pilot1_correlation_records.size() == _pilot2_correlation_records.size())); // ensure both pilot records are at the same point if both are present

    _start_pilot_correlation_records_index = (_pilot1_correlation_records.size() / (_analysis_window_size*_prn_per_symbol)) * _analysis_window_size * _prn_per_symbol; // record beginning of new analysis window in the pilot correlation records
    _start_message_correlation_records_index = _message_correlation_records.size(); // record beginning of new analysis window in the message correlation records
    _pilot1_time_shift_occurences.clear();
    _pilot2_time_shift_occurences.clear();
    _histo_basic_time_shift_occurences.clear();
    _validate_count = 0;
    _analysis_index++;
}


//=================================================================================================
void PilotCorrelationAnalyzer::dump_pilot_correlation_records(std::ostringstream& os, bool alternate) const
{
	std::vector<PilotCorrelationRecord>::const_iterator it;
	std::vector<PilotCorrelationRecord>::const_iterator it_end;

	if (alternate && _two_pilots)
	{
		it = _pilot2_correlation_records.begin();
		it_end = _pilot2_correlation_records.end();
	}
	else
	{
		it = _pilot1_correlation_records.begin();
		it_end = _pilot1_correlation_records.end();
	}

	PilotCorrelationRecord::dump_oneline_banner(os);

	for (; it != it_end; ++it)
	{
		it->dump_oneline(os);
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
		it->dump_line(1.0, os);
	}
}


//=================================================================================================
void PilotCorrelationAnalyzer::dump_histo_time_shift_occurences(std::ostringstream& os) const
{
	std::vector<std::vector<PrnTimeShiftRecord> >::const_iterator h_v_it = _histo_time_shift_occurences_vector.begin();
	const std::vector<std::vector<PrnTimeShiftRecord> >::const_iterator h_v_end = _histo_time_shift_occurences_vector.end();

	for (unsigned int i=0; h_v_it != h_v_end; ++h_v_it, i++)
	{
	    std::vector<PrnTimeShiftRecord>::const_iterator histo_it = h_v_it->begin();
	    const std::vector<PrnTimeShiftRecord>::const_iterator histo_end = h_v_it->end();

	    os << "#" << i << ": Pilot " << _pilot_numbers[i] << ": [";

	    for (; histo_it != histo_end; ++histo_it)
	    {
	        if (histo_it != h_v_it->begin())
	        {
	            os << ", ";
	        }

	        os << "(";
	        histo_it->dump(os, _fft_N);
	        os << ")";
	    }

	    os << "]" << std::endl;
	}
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


//=================================================================================================
PilotCorrelationAnalyzer::PrnTimeShiftRecord::PrnTimeShiftRecord() :
    peak_start(0),
    peak_length(0),
    peak_mag_sum(0),
    peak_shift_at_max(0),
    peak_max(0)
{}


//=================================================================================================
void PilotCorrelationAnalyzer::PrnTimeShiftRecord::dump(std::ostringstream& os, unsigned int time_shifts_size) const
{
    WsgcUtils::print_interval(os, peak_start, peak_length, time_shifts_size);
    os << " :" << peak_shift_at_max << "(" << peak_max << ") : " << peak_mag_sum;
}


//=================================================================================================
bool PilotCorrelationAnalyzer::PrnTimeShiftRecord::order_peak_mag_sum(const PrnTimeShiftRecord& left, const PrnTimeShiftRecord& right)
{
    return left.peak_mag_sum > right.peak_mag_sum;
}
