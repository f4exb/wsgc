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

     Final correlation data analysis and message estimation 
*/
#include "WsgcUtils.h"
#include "DecisionBox_MFSK.h"
#include "MFSK_MessageDemodulationRecord.h"
#include "DecisionBox_Thresholds.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cassert>



//=================================================================================================
DecisionBox_MFSK::DecisionBox_MFSK(
		unsigned int fft_size,
		unsigned int nb_fft_per_symbol,
        const DecisionBox_Thresholds& decision_thresholds) :
	_fft_size(_fft_size),
	_nb_fft_per_symbol(nb_fft_per_symbol),
    _decision_thresholds(decision_thresholds)
{}


//=================================================================================================
DecisionBox_MFSK::~DecisionBox_MFSK()
{}


//=================================================================================================
void DecisionBox_MFSK::estimate_symbols(const std::vector<MFSK_MessageDemodulationRecord>& demodulation_records)
{
	std::vector<MFSK_MessageDemodulationRecord>::const_iterator demod_it = demodulation_records.begin();
	const std::vector<MFSK_MessageDemodulationRecord>::const_iterator demod_end = demodulation_records.end();

	for (; demod_it != demod_end; ++demod_it)
	{
		if ((demod_it->_max_magnitude / demod_it->_avg_magnitude) > _decision_thresholds._peak_margin_threshold)
		{
			_decoded_symbols.push_back(demod_it->_symbol_ordinal);
		}
		else
		{
			_decoded_symbols.push_back(-1);
		}
	}
}

//=================================================================================================
void DecisionBox_MFSK::dump_decision_status(
		std::ostringstream& os,
		std::vector<unsigned int>& original_symbols,
		const std::vector<MFSK_MessageDemodulationRecord>& demodulation_records,
		bool no_trivial) const
{
	assert(_decoded_symbols.size() == demodulation_records.size());
	std::vector<decision_status_t> symbols_decision_status;
    unsigned int count_false_reject = 0;
    unsigned int count_false_accept = 0;
    unsigned int count_true_reject = 0;
    unsigned int count_true_accept = 0;

	for (unsigned int i=0; i < original_symbols.size(); i++)
	{
		if (i < _decoded_symbols.size())
		{
			if (_decoded_symbols[i] == -1)
			{
				if (demodulation_records[i]._symbol_ordinal == original_symbols[i])
				{
					symbols_decision_status.push_back(decision_status_false_reject);
                    count_false_reject++;
				}
				else
				{
					symbols_decision_status.push_back(decision_status_true_reject);
                    count_true_reject++;
				}
			}
			else
			{
				if (_decoded_symbols[i] == original_symbols[i])
				{
					symbols_decision_status.push_back(decision_status_true_accept);
                    count_true_accept++;
				}
				else
				{
					symbols_decision_status.push_back(decision_status_false_accept);
                    count_false_accept++;
				}
			}
		}
	}

	std::vector<decision_status_t>::const_iterator sit = symbols_decision_status.begin();
	const std::vector<decision_status_t>::const_iterator send = symbols_decision_status.end();
	unsigned int si=0;

	for (; sit != send; ++sit, si++)
	{
		if ((*sit != decision_status_true_accept) || !no_trivial)
		{
			os << "_DCS " << std::setw(3) << si << " ";
			dump_decoding_status(os, *sit);
			os << " " << demodulation_records[si]._max_magnitude/demodulation_records[si]._avg_magnitude << std::endl;
		}
	}
    
    unsigned int count_erasures = count_true_reject + count_false_reject;
    unsigned int count_errors = count_false_accept;
    unsigned int rs_corr = count_erasures + 2*count_errors;
    os << "_SDC " << count_true_accept << "," << count_true_reject << "," << count_false_accept << "," << count_false_reject << "," << rs_corr << std::endl;
}


//=================================================================================================
void DecisionBox_MFSK::dump_decoding_status(std::ostringstream& os, decision_status_t decision_status) const
{
	switch (decision_status)
	{
	case decision_status_false_accept:
		os << "FA";
		break;
	case decision_status_false_reject:
		os << "FR";
		break;
	case decision_status_true_accept:
		os << "TA";
		break;
	case decision_status_true_reject:
		os << "TR";
		break;
	default:
		os << "XX";
		break;
	}
}

