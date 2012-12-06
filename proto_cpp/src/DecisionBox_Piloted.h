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
#ifndef __DECISION_BOX_PILOTED_H__
#define __DECISION_BOX_PILOTED_H__

#include "WsgcTypes.h"
#include "DecisionBox.h"
#include <vector>

class PilotCorrelationAnalyzer;

/**
 * \brief Correlation data analysis and message estimation. For operation with pilot(s) PRN(s)
 * transmitted simultaneously with the message PRNs
 */
class DecisionBox_Piloted : public DecisionBox
{
public:
	/**
	 * Constructor
	 * \param prn_per_symbol Number of PRNs per symbol
     * \param fft_N Size of the FFT, this is also the number of samples in one PRN
     * \param pilot_correlation_analyzer Reference to the pilot correlation analyzer
	 */
	DecisionBox_Piloted(unsigned int prn_per_symbol, unsigned int fft_size, const PilotCorrelationAnalyzer& pilot_correlation_analyzer);

    virtual ~DecisionBox_Piloted();

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
        
	static const wsgc_float max_to_avg_threshold;   //<! max / average ratio threshold (for piloted operation)
	static const wsgc_float signal_to_noise_avg_threshold; //<! signal / noise average threshold (for piloted operation)

	bool challenge_matching_symbol(std::vector<CorrelationRecord>::const_iterator& matching_record_it);

    /**
     * Select PRN for the estimation process depending on corresponding correlation record values and the preferred PRN time shift
     * \param record_it Reference to the iterator in the message correlation records vector
     * \return True if selected else false
	*/
	bool select_record(std::vector<CorrelationRecord>::const_iterator& record_it);
};

#endif // __DECISION_BOX_PILOTED_H__
