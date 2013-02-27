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
     In this case the proper synchronization of symbols is assumed and assumed starting
     at a definite PRN index.
*/
#ifndef __DECISION_BOX_PILOTED_AND_SYNCED_H__
#define __DECISION_BOX_PILOTED_AND_SYNCED_H__

#include "WsgcTypes.h"
#include "DecisionBox.h"
#include <vector>

class PilotCorrelationAnalyzer;
class DecisionBox_Thresholds;

/**
 * \brief Correlation data analysis and message estimation. For operation with pilot(s) PRN(s)
 * transmitted simultaneously with the message PRNs and symbol boundary synchronization. It is expected
 * that the symbol starts at a definite PRN index.
 */
class DecisionBox_Piloted_And_Synced : public DecisionBox
{
public:
	/**
	 * Constructor
	 * \param prn_per_symbol Number of PRNs per symbol
     * \param fft_N Size of the FFT, this is also the number of samples in one PRN
     * \param decision_thresholds Reference to the decision thresholds
     * \param pilot_correlation_analyzer Reference to the pilot correlation analyzer
	 */
	DecisionBox_Piloted_And_Synced(unsigned int prn_per_symbol, 
        unsigned int fft_size, 
        const DecisionBox_Thresholds& decision_thresholds,
        const PilotCorrelationAnalyzer& pilot_correlation_analyzer);

    virtual ~DecisionBox_Piloted_And_Synced();

    /**
	 * Does the preparatory analysis phase
	 */
    virtual void analyze_records();

	/**
	 * Does the message symbols estimation
	 */
    virtual void estimate_symbols();
        
protected:
    const PilotCorrelationAnalyzer& _pilot_correlation_analyzer;
        
	bool challenge_matching_symbol(std::vector<CorrelationRecord>::const_iterator& matching_record_it);
	bool test_maxavg_override(std::vector<CorrelationRecord>::const_iterator& matching_record_it, wsgc_float ratio);

    /**
     * Select PRN for the estimation process depending on corresponding correlation record values and the preferred PRN time shift
     * \param record_it Reference to the iterator in the message correlation records vector
     * \return True if selected else false
	*/
	bool select_record(std::vector<CorrelationRecord>::const_iterator& record_it);
};

#endif // __DECISION_BOX_PILOTED_AND_SYNCED_H__
