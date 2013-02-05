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
      
     Decision box

     Final correlation data analysis and message estimation for operation without
     pilot(s) PRN(s) transmitted simultaneously with the message PRNs
     In this case the proper synchronization of symbols is assumed and assumed starting
     at a definite PRN index. There are however only one correlation record per symbol
     sampled at last PRN in symbol cycle.
*/
#include "DecisionBox_Unpiloted_And_Synced.h"
#include "CorrelationRecord.h"

#include <iostream>
#include <iomanip>


//=================================================================================================
DecisionBox_Unpiloted_And_Synced::DecisionBox_Unpiloted_And_Synced(unsigned int prn_per_symbol,
		unsigned int fft_size,
		std::vector<CorrelationRecord>& correlation_records) :
	DecisionBox(prn_per_symbol, fft_size),
	_correlation_records(correlation_records),
	_max_to_avg_ok_threshold(1.1),
	_max_to_avg_ok_threshold_cuda(1.1)
{
	_prni_at_max_invalid = false; // N/A in this case
}


//=================================================================================================
DecisionBox_Unpiloted_And_Synced::~DecisionBox_Unpiloted_And_Synced()
{}


//=================================================================================================
void DecisionBox_Unpiloted_And_Synced::analyze_records()
{
	// Useless in this case
}


//=================================================================================================
void DecisionBox_Unpiloted_And_Synced::estimate_symbols()
{
	std::vector<CorrelationRecord>::const_iterator record_it = _correlation_records.begin();
	const std::vector<CorrelationRecord>::const_iterator record_end = _correlation_records.end();
	unsigned int record_nb = 0;

	for (;record_it != record_end; ++record_it, record_nb++)
	{
        // print current correlation record
        std::ostringstream os;
        os << std::setw(3) << std::setfill('0') << record_nb << " ";
        record_it->dump(_mag_display_adj_factor, os);
        std::cout << os.str() << std::endl;

        DecisionRecord& decision_record = new_decision_record();
        decision_record.symbol_index = record_nb; // exactly one record per symbol
        decision_record.global_prn_index = record_it->global_prn_index;
        decision_record.prn_per_symbol_index = record_it->global_prn_index % _prn_per_symbol;
        decision_record.prn_index_max = record_it->prn_index_max;
        decision_record.select_count = 1;
        decision_record.magnitude_max = record_it->magnitude_max;
        decision_record.magnitude_avg = record_it->magnitude_avg;
        decision_record.noise_avg = record_it->noise_avg;
        decision_record.shift_index_max = record_it->shift_index_max;
        decision_record.f_rx = record_it->f_rx;

		if (record_it->selected)
		{
			if (challenge_matching_symbol(record_it))
			{
				decision_record.decision_type = DecisionRecord::decision_ok_strong;
				decision_record.validated = true;
				_decoded_symbols.push_back(record_it->prn_index_max);
			}
			else
			{
				decision_record.decision_type = DecisionRecord::decision_ko_too_weak;
				decision_record.validated = false;
				_decoded_symbols.push_back(-1); // not valid => erasure symbol
			}
		}
		else
		{
			decision_record.decision_type = DecisionRecord::decision_ko_no_valid_rec; // by extension
			decision_record.validated = false;
			_decoded_symbols.push_back(-1); // not valid => erasure symbol
		}
	}
}
        

//=================================================================================================
bool DecisionBox_Unpiloted_And_Synced::challenge_matching_symbol(std::vector<CorrelationRecord>::const_iterator& matching_record_it)
{
	if (matching_record_it->magnitude_avg == 0)
	{
		return true; // disabled
	}

    wsgc_float max_to_avg = matching_record_it->magnitude_max / matching_record_it->magnitude_avg;
	std::cout << "    +  challenge rec # " << matching_record_it->global_prn_index << std::endl;

    if (max_to_avg > (_use_cuda ? _max_to_avg_ok_threshold_cuda : _max_to_avg_ok_threshold))
    {
    	return true;
    }
    else
    {
    	return false;
    }
}


//=================================================================================================
void DecisionBox_Unpiloted_And_Synced::set_thresholds(Modulation& modulation)
{
	if (modulation.getScheme() == Modulation::Modulation_OOK)
	{
		_max_to_avg_ok_threshold_cuda = 1.06;
		_max_to_avg_ok_threshold = 1.06;
	}
	else if (modulation.getScheme() == Modulation::Modulation_DBPSK)
	{
		_max_to_avg_ok_threshold_cuda = 1.15;
		_max_to_avg_ok_threshold = 1.15;
	}
	// else leave defaults (normally not used)
}
