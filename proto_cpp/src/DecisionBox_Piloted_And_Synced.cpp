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
#include "DecisionBox_Piloted_And_Synced.h"
#include "PilotCorrelationAnalyzer.h"
#include "DecisionRecord.h"
#include "DecisionBox_Thresholds.h"
#include <iostream>
#include <iomanip>
#include <algorithm>


//=================================================================================================
DecisionBox_Piloted_And_Synced::DecisionBox_Piloted_And_Synced(unsigned int prn_per_symbol, 
        unsigned int fft_size, 
        const DecisionBox_Thresholds& decision_thresholds,
        const PilotCorrelationAnalyzer& pilot_correlation_analyzer) :
    DecisionBox(prn_per_symbol, fft_size, decision_thresholds),
    _pilot_correlation_analyzer(pilot_correlation_analyzer)
{
	std::cout << "using DecisionBox_Piloted_And_Synced" << std::endl;
    _preferred_symbol_prn_i = _prn_per_symbol-1;
    _prni_at_max_invalid = false;

}


//=================================================================================================
DecisionBox_Piloted_And_Synced::~DecisionBox_Piloted_And_Synced()
{}


//=================================================================================================
void DecisionBox_Piloted_And_Synced::analyze_records()
{
    // Useless in this case:
	//DecisionBox::estimate_preferred_symbol_prn_i(_pilot_correlation_analyzer.get_message_correlation_records()); // find preferred PRN per symbol index (symbol synchronization)
}


//=================================================================================================
void DecisionBox_Piloted_And_Synced::estimate_symbols()
{
    unsigned int nb_symbols = _pilot_correlation_analyzer.get_message_correlation_records().size() / _prn_per_symbol;
    //unsigned int preferred_symbol_prn_i;
    //unsigned int preferred_shift;
    
    if (nb_symbols < 4)
    {
        std::cout << "DecisionBox: not enough symbols to do a proper estimation" << std::endl;
    }
    else
    {
        std::vector<CorrelationRecord>::const_iterator record_it = _pilot_correlation_analyzer.get_message_correlation_records().begin();
        const std::vector<CorrelationRecord>::const_iterator records_end = _pilot_correlation_analyzer.get_message_correlation_records().end();
        std::vector<CorrelationRecord>::const_iterator matching_record_it;
        
        unsigned int record_i = 0; // always start at 0 for synced message
        unsigned int symbol_cycle_nb;
        unsigned int prn_index_in_cycle;
        unsigned int select_count;
        
        unsigned int record_nb = 0;
        bool matching_at_preferred_index;
        
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
                matching_record_it = records_end; // reset matching record iterator (not found)
                select_count = 0;
            }
            
			// try to find eligible record
			if (select_record(record_it))
			{
	            // do not store first in cycle (low confidence)
	            if (prn_index_in_cycle > 0)
	            {
					matching_record_it = record_it;
					matching_at_preferred_index = (record_it->prn_per_symbol_index == _preferred_symbol_prn_i); // preferred PRN index in symbol (maximum integration)
	            }

				select_count++;
            }
            
            // take decision on last PRN of the current symbol cycle
            if (prn_index_in_cycle == _prn_per_symbol - 1) 
            {
                DecisionRecord& decision_record = new_decision_record();
            
                if (matching_record_it == records_end)
                {
                    std::cout << "    +  no valid PRN found in cycle" << std::endl;
                    decision_record.decision_type = DecisionRecord::decision_ko_no_valid_rec;
                    _decoded_symbols.push_back(-1); // nothing found => erasure symbol 
                }
                else 
                {
                    decision_record.symbol_index = matching_record_it->global_prn_index / _prn_per_symbol;
                    decision_record.global_prn_index = matching_record_it->global_prn_index;
                    decision_record.prn_per_symbol_index = matching_record_it->global_prn_index % _prn_per_symbol;
                    decision_record.prn_index_max = matching_record_it->prn_index_max;
                    decision_record.select_count = select_count;
                    decision_record.magnitude_max = matching_record_it->magnitude_max;
                    decision_record.magnitude_avg = matching_record_it->magnitude_avg;
                    decision_record.noise_avg = matching_record_it->noise_avg;
                    decision_record.shift_index_max = matching_record_it->pilot_shift;
                    decision_record.f_rx = matching_record_it->f_rx;
                    
                    if (select_count < (3 * _prn_per_symbol) / 4)
                    {
                    	if (test_maxavg_override(matching_record_it, (wsgc_float)select_count / (wsgc_float)_prn_per_symbol))
                    	{
							std::cout << "    +  symbol validated at end of cycle even if there are not enough selected records for this cycle" << std::endl;
							decision_record.decision_type = DecisionRecord::decision_ok_not_enough_rec;
							decision_record.validated = true;
							_decoded_symbols.push_back(matching_record_it->prn_index_max);
                    	}
                    	else
                    	{
							std::cout << "    +  symbol rejected at end of cycle because there are not enough selected records for this cycle" << std::endl;
							decision_record.decision_type = DecisionRecord::decision_ko_not_enough_rec;
							_decoded_symbols.push_back(-1); // not valid => erasure symbol
                    	}
                    }
                    else if (challenge_matching_symbol(matching_record_it))
                    {
                        if (matching_at_preferred_index)
                        {
                            std::cout << "    +  symbol validated at preferred index" << std::endl;
                        }
                        else
                        {
                            std::cout << "    +  symbol validated not at preferred index" << std::endl;
                        }
                        
                        decision_record.decision_type = DecisionRecord::decision_ok_strong;
                        decision_record.validated = true;
                        _decoded_symbols.push_back(matching_record_it->prn_index_max);
                    }
                    else
                    {
                        std::cout << "    +  symbol not validated because it is too weak" << std::endl;
                        decision_record.decision_type = DecisionRecord::decision_ko_too_weak;
                        _decoded_symbols.push_back(-1); // not valid => erasure symbol 
                    }
                }
            }
            
            /*
                // preferred shift encountered
                // challenged last condition because it skips the next identical symbol if this occurs in a sequence. 
                // it can be validated at preferred index and skipped at other indexes
                // if ((record_it->shift_index_max == preferred_shift) && (last_validated_symbol != record_it->prn_index_max)) 
                //if ((record_it->shift_index_max == preferred_shift)
                //   && ((record_it->magnitude_max / record_it->noise_max) > noise_margin_threshold)) // and if maximum peak / maximum noise is above threshold
                if (select_record(record_it))
                {
                    //if (symbol_found_for_cycle_nb != symbol_cycle_nb) // do not process if a valid symbol has already been found in this cycle
                    //{
                        std::ostringstream os;
                        os << std::setiosflags(std::ios_base::fixed);
                        
                        if (record_it->prn_per_symbol_index == _preferred_symbol_prn_i) // and preferred PRN index in symbol (maximum confidence)
                        {
                        	if (select_count < _prn_per_symbol / 2)
                        	//if ((prn_index_in_cycle == _prn_per_symbol - 1) && (select_count < _prn_per_symbol / 2))
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
                    //}
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
						//if (symbol_matching_not_at_preferred_index != last_validated_symbol)
						//{
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
						//}
                	}
                }
                
                last_prn_index_max = record_it->prn_index_max;
            }
            */
            
            record_i++;
        }
    }
}


//=================================================================================================
bool DecisionBox_Piloted_And_Synced::select_record(std::vector<CorrelationRecord>::const_iterator& record_it)
{
	return (record_it->selected);
}


//=================================================================================================
bool DecisionBox_Piloted_And_Synced::challenge_matching_symbol(std::vector<CorrelationRecord>::const_iterator& matching_record_it)
{
	if (matching_record_it->magnitude_avg == 0)
	{
		return true; // disabled
	}

    wsgc_float max_to_avg = matching_record_it->magnitude_max / matching_record_it->magnitude_avg;
	std::cout << "    +  challenge rec # " << matching_record_it->global_prn_index << std::endl;
    
    if (matching_record_it->noise_avg == 0)
    {
        return max_to_avg > _decision_thresholds._max_to_avg_ok_threshold; // noise test disabled
    }
    else
    {
        wsgc_float signal_to_noise = matching_record_it->magnitude_max / matching_record_it->noise_avg;
        if (signal_to_noise < _decision_thresholds._signal_to_noise_avg_ko_threshold)
        {
        	return false;
        }
        else if (max_to_avg > _decision_thresholds._max_to_avg_ok_threshold)
        {
        	return true;
        }
        else if (max_to_avg > _decision_thresholds._max_to_avg_cdt_threshold)
        {
        	if (signal_to_noise > _decision_thresholds._signal_to_noise_avg_cdt_threshold)
        	{
        		return true;
        	}
        	else
        	{
        		return false;
        	}
        }
        else
        {
        	return false;
        }
    }
}


//=================================================================================================
bool DecisionBox_Piloted_And_Synced::test_maxavg_override(
		std::vector<CorrelationRecord>::const_iterator& matching_record_it,
		wsgc_float ratio)
{
	if (matching_record_it->magnitude_avg == 0)
	{
		return true; // disabled
	}

    wsgc_float max_to_avg = matching_record_it->magnitude_max / matching_record_it->magnitude_avg;
	std::cout << "    +  maxavg override test rec # " << matching_record_it->global_prn_index << std::endl;

	if (ratio < 0.5)
	{
		return max_to_avg > _decision_thresholds._max_to_avg_ok_threshold_ner_0_5;
	}
	else if (ratio < 0.75)
	{
		return max_to_avg > _decision_thresholds._max_to_avg_ok_threshold_ner_0_75;
	}
	else
	{
		return true;
	}
}

