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
#include "DecisionBox_Unpiloted.h"
#include "Modulation.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>

const wsgc_float DecisionBox_Unpiloted::noise_margin_threshold = 1.3; // max / noise max ratio threshold. Used to be 1.2
const wsgc_float DecisionBox_Unpiloted::max_to_avg_margin_threshold = 1.1; // max / avg magnitude ratio threshold.

//=================================================================================================
DecisionBox_Unpiloted::DecisionBox_Unpiloted(unsigned int prn_per_symbol,
		unsigned int fft_size,
		const std::vector<CorrelationRecord>& correlation_records,
		const std::map<unsigned int, unsigned int>& shift_occurences,
		const Modulation& modulation) :
    DecisionBox(prn_per_symbol, fft_size),
    _correlation_records(correlation_records),
    _shift_occurences(shift_occurences),
    _modulation(modulation)
{}

  
//=================================================================================================
DecisionBox_Unpiloted::~DecisionBox_Unpiloted()
{}


//=================================================================================================
void DecisionBox_Unpiloted::analyze_records()
{
    std::map<unsigned int, unsigned int>::const_iterator shift_occurences_it;
    
    DecisionBox::estimate_preferred_symbol_prn_i(_correlation_records); // find preferred PRN per symbol index (symbol synchronization)
    
    if (_correlation_records.size() > std::max(_prn_per_symbol,2u))
    {
        for (shift_occurences_it = _shift_occurences.begin(); shift_occurences_it != _shift_occurences.end(); ++shift_occurences_it)
        {
            _histo_shift_occurences.push_back(*shift_occurences_it);
        }
        
        std::sort(_histo_shift_occurences.begin(), _histo_shift_occurences.end(), DecisionBox::histo_order);
        
        std::ostringstream shift_os;
        WsgcUtils::print_histo(_histo_shift_occurences, shift_os);
        
        std::cout << "Shifts .............................: " << shift_os.str() << std::endl;
    }
    else
    {
        std::cout << "DecisionBox: not enough correlation records to analyze best PRN time shift" << std::endl;
    }
}


//=================================================================================================
void DecisionBox_Unpiloted::estimate_symbols()
{
    unsigned int nb_symbols = _correlation_records.size() / _prn_per_symbol;
    unsigned int preferred_symbol_prn_i;
    unsigned int preferred_shift;
    
    if (nb_symbols < 4)
    {
        std::cout << "DecisionBox: not enough symbols to do a proper estimation" << std::endl;
    }
    else if (_shift_occurences.size() == 0)
    {
        std::cout << "DecisionBox: no shift occurence history" << std::endl;
    }
    else if (_symbol_prn_i_at_max.size() == 0)
    {
        std::cout << "DecisionBox: no PRN in symbol index occurence history" << std::endl;
    }
    else
    {
        //TODO: aggregate neighboring bins using peak dithering data increasing peak histogram weight - done with piloted
        //TODO: among equal candidates select shift that also best matches preferred symbol PRN index
        preferred_shift = _histo_shift_occurences.begin()->first;
        unsigned int preferred_shift_margin = _histo_shift_occurences.begin()->second - (_histo_shift_occurences.begin()+1)->second;
        
        if (((float)preferred_shift_margin/_correlation_records.size()) < peak_margin_threshold)
        {
            std::cout << "Best maximum peak margin is too small: message should be invalidated" << std::endl;
        }
        
        std::cout << "Preferred shift is " << preferred_shift << std::endl;
        
        std::vector<CorrelationRecord>::const_iterator record_it = _correlation_records.begin();
        const std::vector<CorrelationRecord>::const_iterator records_end = _correlation_records.end();
        std::vector<CorrelationRecord>::const_iterator matching_record_it;
        
        unsigned int start_of_cycle_index = (_preferred_symbol_prn_i+1) %  _prn_per_symbol; // make sure cycle ends at //the next index after// the preferred shift, //hence starts two indexes later
        std::cout << "Start of cycle index: " << start_of_cycle_index << std::endl;
        unsigned int record_i = _prn_per_symbol - start_of_cycle_index; //= 0; position of 0 in the cycle
        unsigned int symbol_cycle_nb;
        unsigned int prn_index_in_cycle;
        
        unsigned int last_prn_index_max = -1;
        unsigned int last_validated_symbol = -2; // init to erasure different from below to prevent possible start condition conflicts
        int symbol_matching_not_at_preferred_index = -1;
        unsigned int symbol_found_for_cycle_nb = -1;
        unsigned int record_nb = 0;
        bool started = false;
        
        for (; record_it != records_end; ++record_it, record_nb++)
        {
            // print current correlation record
            std::ostringstream os;
            os << std::setw(3) << std::setfill('0') << record_nb / _prn_per_symbol << " ";
            //record_it->dump(_fft_size * _prn_per_symbol, os);
            record_it->dump(_mag_display_adj_factor, os);
            
            symbol_cycle_nb = record_i / _prn_per_symbol; 
            prn_index_in_cycle = record_i % _prn_per_symbol;
            os << "    -- " << symbol_cycle_nb << ":" << prn_index_in_cycle;
            std::cout << os.str() << std::endl;
            os.flush();
            
            // start of cycle: append erasure symbol and reset symbol found in the cycle outside the preferred index
            if (prn_index_in_cycle == 0)
            {
                if (!started)
                {
                    std::cout << "Synchronize start of cycle 0" << std::endl;
                    record_i = 0;
                    started = true; 
                }
                
                _decoded_symbols.push_back(-1); // erasure symbol by default
                symbol_matching_not_at_preferred_index = -1; // reset symbol
                matching_record_it = record_it;
            }
            
            // process only from the first start and prevent symbol spillover 
            if ((started) && (prn_index_in_cycle > 1))
            {
                // preferred shift encountered
                // challenged last condition because it skips the next identical symbol if this occurs in a sequence. 
                // it can be validated at preferred index and skipped at other indexes
                // if ((record_it->shift_index_max == preferred_shift) && (last_validated_symbol != record_it->prn_index_max)) 
                //if ((record_it->shift_index_max == preferred_shift)
                //   && ((record_it->magnitude_max / record_it->noise_max) > noise_margin_threshold)) // and if maximum peak / maximum noise is above threshold
                if (select_record(record_it, preferred_shift))
                {
                    if (symbol_found_for_cycle_nb != symbol_cycle_nb) // do not process if a valid symbol has already been found in this cycle
                    {
                        std::ostringstream os;
                        os << std::setiosflags(std::ios_base::fixed);
                        
                        if (record_it->prn_per_symbol_index == _preferred_symbol_prn_i) // and preferred PRN index in symbol (maximum confidence)
                        {
                            last_validated_symbol = record_it->prn_index_max;
                            symbol_found_for_cycle_nb = symbol_cycle_nb;
                            
                            _decoded_symbols.back() = record_it->prn_index_max;
                            os << "    +  symbol validated at preferred index" << std::endl;
                        }
                        else // not preferred index (less confidence, defer until possible best match)
                        {
                            if (last_validated_symbol != record_it->prn_index_max) // skip next identical symbol not at preferred index
                            {
                                symbol_matching_not_at_preferred_index = record_it->prn_index_max;
                                matching_record_it = record_it;
                                os << "    +  symbol not matching at preferred index" << std::endl;
                                matching_record_it = record_it; // retain record pointer for end of cycle decision
                            }
                        }
                        
                        std::cout << os.str();
                    }
                    else
                    {
                    	std::cout << "    -  symbol already validated skipping..." << std::endl;
                    }
                }
                else
                {
                	std::cout << "    -  record not selected skipping..." << std::endl;
                }
                
                // end of cycle: if no matching symbol was found take any symbol found in the cycle outside the preferred index (that would othewise have matched)
                if ((prn_index_in_cycle == _prn_per_symbol - 1) && (symbol_matching_not_at_preferred_index > 0) && (symbol_found_for_cycle_nb != symbol_cycle_nb))
                {
                    if (symbol_matching_not_at_preferred_index != last_validated_symbol)
                    {
                        last_validated_symbol = symbol_matching_not_at_preferred_index;
                        
                        _decoded_symbols.back() = symbol_matching_not_at_preferred_index;
                        std::cout << "    +  symbol validated at end of cycle " << symbol_cycle_nb << std::endl;
                    }
                }
                
                last_prn_index_max = record_it->prn_index_max;
            }
            
            record_i++;
        }
    }
}


//=================================================================================================
bool DecisionBox_Unpiloted::select_record(std::vector<CorrelationRecord>::const_iterator& record_it, unsigned int preferred_shift)
{
	if (_modulation.getScheme() == Modulation::Modulation_OOK)
	{
		return ((record_it->shift_index_max == preferred_shift) && ((record_it->magnitude_max / record_it->magnitude_avg) > max_to_avg_margin_threshold));
	}
	else if (_modulation.getScheme() == Modulation::Modulation_BPSK)
	{
		return ((record_it->shift_index_max == preferred_shift) && ((record_it->magnitude_max / record_it->noise_max) > noise_margin_threshold));
	}
	else
	{
		return record_it->shift_index_max == preferred_shift;
	}
}
