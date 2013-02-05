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

     Time correlation analyzer

     Analyzes time shift data of correlation records.
*/

#include "TimeCorrelationAnalyzer.h"
#include "WsgcUtils.h"
#include <algorithm>

const wsgc_float TimeCorrelationAnalyzer_Common::_best_time_margin_threshold = 0.13; // selected correlation time delta difference with next / number of correlations ratio threshold

//=================================================================================================
TimeCorrelationAnalyzer_Common::TimeCorrelationAnalyzer_Common(unsigned int fft_N, unsigned int time_index_tolerance) :
    _fft_N(fft_N),
    _time_index_tolerance(time_index_tolerance), //TODO: use time index tolerance
    _time_correlations_count(0),
    _sum_corr_max(0.0)
{
	std::cout << "TimeCorrelationAnalyzer_Common::TimeCorrelationAnalyzer_Common" << std::endl;
}


//=================================================================================================
TimeCorrelationAnalyzer_Common::~TimeCorrelationAnalyzer_Common()
{
}


//=================================================================================================
void TimeCorrelationAnalyzer_Common::reset_analysis()
{
    _time_shift_occurences.clear();
    _histo_basic_time_shift_occurences.clear();
    _time_correlations_count = 0;
    _sum_corr_max = 0.0;
}


//=================================================================================================
void TimeCorrelationAnalyzer_Common::dump_histo_time_shift_occurences(std::ostringstream& os) const
{
	std::vector<std::vector<PrnTimeShiftRecord> >::const_iterator h_v_it = _histo_time_shift_occurences_vector.begin();
	const std::vector<std::vector<PrnTimeShiftRecord> >::const_iterator h_v_end = _histo_time_shift_occurences_vector.end();

	for (unsigned int i=0; h_v_it != h_v_end; ++h_v_it, i++)
	{
	    std::vector<PrnTimeShiftRecord>::const_iterator histo_it = h_v_it->begin();
	    const std::vector<PrnTimeShiftRecord>::const_iterator histo_end = h_v_it->end();

	    os << "#" << i << ": [";

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
void TimeCorrelationAnalyzer_Common::make_time_shift_histogram()
{
	static const std::vector<PrnTimeShiftRecord> tmp_histo_time_shift_occurences;
    _histo_time_shift_occurences_vector.push_back(tmp_histo_time_shift_occurences);

    std::map<unsigned int, unsigned int>::const_iterator shift_occurences_it = _time_shift_occurences.begin();
    const std::map<unsigned int, unsigned int>::const_iterator shift_occurences_end = _time_shift_occurences.end();
    std::map<unsigned int, unsigned int>::const_iterator last_shift_occurence_it = _time_shift_occurences.begin();
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
TimeCorrelationAnalyzer_Common::PrnTimeShiftRecord::PrnTimeShiftRecord() :
    peak_start(0),
    peak_length(0),
    peak_mag_sum(0),
    peak_shift_at_max(0),
    peak_max(0),
    nb_corr_records(0)
{}


//=================================================================================================
void TimeCorrelationAnalyzer_Common::PrnTimeShiftRecord::dump(std::ostringstream& os, unsigned int time_shifts_size) const
{
    WsgcUtils::print_interval(os, peak_start, peak_length, time_shifts_size);
    os << " :" << peak_shift_at_max << "(" << peak_max << ") : " << peak_mag_sum;
}


//=================================================================================================
bool TimeCorrelationAnalyzer_Common::PrnTimeShiftRecord::order_peak_mag_sum(const PrnTimeShiftRecord& left, const PrnTimeShiftRecord& right)
{
    return left.peak_mag_sum > right.peak_mag_sum;
}

