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
#ifndef __DECISION_BOX_UNPILOTED_AND_SYNCED_H__
#define __DECISION_BOX_UNPILOTED_AND_SYNCED_H__

#include "WsgcTypes.h"
#include "DecisionBox.h"
#include "Modulation.h"
#include <vector>

class CorrelationRecord;

/**
 * \brief Correlation data analysis and message estimation.
 * For operation without pilot(s) PRN(s) transmitted simultaneously with the message PRNs
 * In this case the proper synchronization of symbols is assumed and assumed starting at
 * a definite PRN index. There are however only one correlation record per symbol sampled
 * at last PRN in symbol cycle.
 */
class DecisionBox_Unpiloted_And_Synced : public DecisionBox
{
public:
	/**
	 * Constructor
	 * \param prn_per_symbol Number of PRNs per symbol
     * \param fft_N Size of the FFT, this is also the number of samples in one PRN
     * \param correlation_records Reference to the vector of correlation records
	 */
	DecisionBox_Unpiloted_And_Synced(unsigned int prn_per_symbol,
			unsigned int fft_size,
			std::vector<CorrelationRecord>& correlation_records);

    virtual ~DecisionBox_Unpiloted_And_Synced();

    /**
	 * Does the preparatory analysis phase
	 */
    virtual void analyze_records();

	/**
	 * Does the message symbols estimation
	 */
    virtual void estimate_symbols();

    /**
     * Sets thresholds depending on modulation type
     */
    void set_thresholds(Modulation& modulation);
        
protected:
	wsgc_float _max_to_avg_ok_threshold;   //<! max / average ratio threshold of unconditionnal acceptance
	wsgc_float _max_to_avg_ok_threshold_cuda; //<! max / average ratio threshold of unconditionnal acceptance for CUDA version
    std::vector<CorrelationRecord>& _correlation_records;

	bool challenge_matching_symbol(std::vector<CorrelationRecord>::const_iterator& matching_record_it);
};

#endif // __DECISION_BOX_UNPILOTED_AND_SYNCED_H__
