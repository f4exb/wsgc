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

     Final correlation data analysis and message estimation.
     For operation when only the message PRNs are transmitted (no pilot)
*/
#ifndef __DECISION_BOX_UNPILOTED_H__
#define __DECISION_BOX_UNPILOTED_H__

#include "DecisionBox.h"
#include "WsgcTypes.h"
#include "CorrelationRecord.h"
#include <vector>
#include <map>

/**
 * \brief SampleAnalyzer Class to handle PRN samples sequencing for operation with only the message PRNs
 */
class DecisionBox_Unpiloted : public DecisionBox
{
    public:
		/**
		 * Constructor
		 * \param prn_per_symbol Number of PRNs per symbol
         * \param fft_N Size of the FFT, this is also the number of samples in one PRN
         * \param correlation_records Reference to the correlation records obtained from the message correlator
         * \param shift_occurences References to the PRN shift occurences 
		 */
        DecisionBox_Unpiloted(unsigned int prn_per_symbol, unsigned int fft_size, const std::vector<CorrelationRecord>& correlation_records, const std::map<unsigned int, unsigned int>& shift_occurences);
        
        virtual ~DecisionBox_Unpiloted();

		/**
		 * Does the preparatory analysis phase
		 */
        virtual void analyze_records();
        
		/**
		 * Does the message symbols estimation 
		 */
        virtual void estimate_symbols();
        
    protected:
        const std::vector<CorrelationRecord>& _correlation_records; //!< Reference to the calculated correlation records
        const std::map<unsigned int, unsigned int>& _shift_occurences; //!< Reference to the calculated PRN time shifts occurences
        std::vector<std::pair<unsigned int, unsigned int> > _histo_shift_occurences; //!< Histogram of the PRN time shifts occurences

        static const wsgc_float noise_margin_threshold; //<! max / noise max ratio threshold

        /**
         * Select PRN for the estimation process depending on corresponding correlation record values and the preferred PRN time shift
         * \param record_it Reference to the iterator in the message correlation records vector
         * \param preferred_shift Estimated preferred PRN time shift
         * \return True if selected else false
    	*/
        bool select_record(std::vector<CorrelationRecord>::const_iterator& record_it, unsigned int preferred_shift); //!< Decide if a correlation record should be selected for estimation
};

#endif // __DECISION_BOX_H__
