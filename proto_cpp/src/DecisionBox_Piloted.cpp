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

     Final correlation data analysis and message estimation for operation with
     pilot(s) PRN(s) transmitted simultaneously with the message PRNs
*/
#include "WsgcUtils.h"
#include "DecisionBox_Piloted.h"
#include "PilotCorrelationAnalyzer.h"
#include <iostream>
#include <iomanip>
#include <algorithm>

const wsgc_float DecisionBox_Piloted::max_to_avg_threshold = 2.0; // 3.2
const wsgc_float DecisionBox_Piloted::signal_to_noise_avg_threshold = 2.0; // 1.2

//=================================================================================================
DecisionBox_Piloted::DecisionBox_Piloted(unsigned int prn_per_symbol, unsigned int fft_size, const PilotCorrelationAnalyzer& pilot_correlation_analyzer) :
    DecisionBox(prn_per_symbol, fft_size),
    _pilot_correlation_analyzer(pilot_correlation_analyzer)
{}


//=================================================================================================
DecisionBox_Piloted::~DecisionBox_Piloted()
{}


//=================================================================================================
void DecisionBox_Piloted::analyze_records()
{
	DecisionBox::estimate_preferred_symbol_prn_i(_pilot_correlation_analyzer.get_message_correlation_records()); // find preferred PRN per symbol index (symbol synchronization)
}


//=================================================================================================
void DecisionBox_Piloted::estimate_symbols()
{
    unsigned int nb_symbols = _pilot_correlation_analyzer.get_message_correlation_records().size() / _prn_per_symbol;
    //unsigned int preferred_symbol_prn_i;
    //unsigned int preferred_shift;
    
    if (nb_symbols < 4)
    {
        std::cout << "DecisionBox: not enough symbols to do a proper estimation" << std::endl;
        
    }
    else if (_symbol_prn_i_at_max.size() == 0)
    {
        std::cout << "DecisionBox: no PRN in symbol index occurence history" << std::endl;
    }
    else
    {
        std::vector<CorrelationRecord>::const_iterator record_it = _pilot_correlation_analyzer.get_message_correlation_records().begin();
        const std::vector<CorrelationRecord>::const_iterator records_end = _pilot_correlation_analyzer.get_message_correlation_records().end();
        std::vector<CorrelationRecord>::const_iterator matching_record_it;
        
        unsigned int start_of_cycle_index = (_preferred_symbol_prn_i+1) %  _prn_per_symbol; // make sure cycle ends at //the next index after// the preferred shift, //hence starts two indexes later
        std::cout << "Start of cycle index: " << start_of_cycle_index << std::endl;
        unsigned int record_i = _prn_per_symbol - start_of_cycle_index; //= 0; position of 0 in the cycle
        unsigned int symbol_cycle_nb;
        unsigned int prn_index_in_cycle;
        unsigned int select_count;
        
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
            record_it->dump(_fft_size/2.0, os);
            
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
                select_count = (record_it->selected ? 1 : 0);
            }
            else
            {
            	select_count += (record_it->selected ? 1 : 0);
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
                if (select_record(record_it))
                {
                    if (symbol_found_for_cycle_nb != symbol_cycle_nb) // do not process if a valid symbol has already been found in this cycle
                    {
                        std::ostringstream os;
                        os << std::setiosflags(std::ios_base::fixed);
                        
                        if (record_it->prn_per_symbol_index == _preferred_symbol_prn_i) // and preferred PRN index in symbol (maximum confidence)
                        {
                        	if ((prn_index_in_cycle == _prn_per_symbol - 1) && (select_count < _prn_per_symbol / 2))
                        	{
                        		std::cout << "    +  symbol rejected at preferred index and end of cycle because there are not enough selected records for this cycle" << std::endl;
                        	}
                        	else
                        	{
								last_validated_symbol = record_it->prn_index_max;
								symbol_found_for_cycle_nb = symbol_cycle_nb;

								if (challenge_matching_symbol(record_it))
								{
									_decoded_symbols.back() = record_it->prn_index_max;
									os << "    +  symbol validated at preferred index" << std::endl;
								}
								else
								{
									os << "    +  symbol not validated at preferred index because it is too weak" << std::endl;
								}
                        	}
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
                }
                
                // end of cycle: if no matching symbol was found take any symbol found in the cycle outside the preferred index (that would othewise have matched)
                if ((prn_index_in_cycle == _prn_per_symbol - 1) && (symbol_matching_not_at_preferred_index > 0) && (symbol_found_for_cycle_nb != symbol_cycle_nb))
                {
                	if (select_count < _prn_per_symbol / 2)
                	{
                		std::cout << "    +  symbol rejected at end of cycle because there are not enough selected records for this cycle" << std::endl;
                	}
                	else
                	{
						if (symbol_matching_not_at_preferred_index != last_validated_symbol)
						{
							last_validated_symbol = symbol_matching_not_at_preferred_index;

							if (challenge_matching_symbol(matching_record_it))
							{
								_decoded_symbols.back() = symbol_matching_not_at_preferred_index;
								std::cout << "    +  symbol validated at end of cycle " << symbol_cycle_nb << std::endl;
							}
							else
							{
								std::cout << "    +  symbol rejected at end of cycle " << symbol_cycle_nb << " because it is too weak" << std::endl;
							}
						}
                	}
                }
                
                last_prn_index_max = record_it->prn_index_max;
            }
            
            record_i++;
        }
    }
}


//=================================================================================================
bool DecisionBox_Piloted::select_record(std::vector<CorrelationRecord>::const_iterator& record_it)
{
	return (record_it->selected);
}


//=================================================================================================
bool DecisionBox_Piloted::challenge_matching_symbol(std::vector<CorrelationRecord>::const_iterator& matching_record_it)
{
	std::cout << "    +  challenge rec # " << matching_record_it->global_prn_index << std::endl;

	if (matching_record_it->magnitude_avg == 0)
	{
		return true; // disabled
	}

    wsgc_float max_to_avg = matching_record_it->magnitude_max / matching_record_it->magnitude_avg;
    
    if (matching_record_it->noise_avg == 0)
    {
        return max_to_avg > max_to_avg_threshold; // noise test disabled
    }
    else
    {
        wsgc_float signal_to_noise = matching_record_it->magnitude_max / matching_record_it->noise_avg;
        return (max_to_avg > max_to_avg_threshold) && (signal_to_noise > signal_to_noise_avg_threshold);
    }
}
