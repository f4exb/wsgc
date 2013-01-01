/*
 * PilotedMultiplePrnCorrelator.h
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

     PilotedMultiplePrnCorrelator

     Processes the correlation of multiple possible PRNs using one (or two) pilot PRN(s)
     It works with a pilot correlator and a message correlator referenced from outside.
     - The pilot correlator looks for a correlation on pilot PRN(s) only and produces the information
       on the best possible frequency shift from zero IF and time delta of the PRN origin in the receiver time reference
     - The message correlator uses this information to look for the one among possible PRNs (the PRN "alphabet") that has a distinctive
       peak at the frequency shift and time delta.

 */

#ifndef __PILOTED_MULTIPLE_PRN_CORRELATOR_H__
#define __PILOTED_MULTIPLE_PRN_CORRELATOR_H__

#include "WsgcTypes.h"
#include <vector>
#include <sstream>
#include <time.h>

class PilotCorrelationAnalyzer;
class CorrelationRecord;
class PilotCorrelator;
class PilotedMessageCorrelator;
class PrnAutocorrelator;

/**
 * \brief Correlator engine to get the correlation estimation of PRNs sent in the message using a pilot correlator and a message correlator
 *
 *  It works with a pilot correlator and a message correlator referenced from outside.
 *  - The pilot correlator looks for a correlation on pilot PRN(s) only and produces the information
 *    on the best possible frequency shift from zero IF and time delta of the PRN origin in the receiver time reference
 *  - The message correlator uses this information to look for the one among possible PRNs (the PRN "alphabet") that has a distinctive
 *    peak at the frequency shift and time delta.
 *
 *  Supports pipelined processing in order to better accomodate the CUDA implementation.
 *  - Pipeline length is set by the Pilot Correlator "PRN batch factor" which can be no less than the averaging length (in PRN numbers) minus one. The
 *    pipeline length is calculated by the Pilot Correlator and is available through a getter method
 *  - is_correlation_record_available method tells if a correlation record is available after execution
 *  - is_correlation_record_available_next method tells if a correlation record will be available next. i.e. is false if it is the last correlation record.
 *
 */
class PilotedMultiplePrnCorrelator
{
public:
    
    /**
    * Correlator engine constructor
    * \param pilot_correlation_analyzer Reference to the pilot correlation analyzer
    * \param correlation_records Reference to the (message) correlation records
    * \param pilot_correlator Reference to the pilot correlator
    * \param message_correlator Reference to the message correlator
    * \param prn_autocorrelator Reference to the PRN autocorrelator
    */
    PilotedMultiplePrnCorrelator(
            PilotCorrelationAnalyzer& pilot_correlation_analyzer,
    		std::vector<CorrelationRecord>& correlation_records,
            PilotCorrelator& pilot_correlator, 
            PilotedMessageCorrelator& message_correlator,
            PrnAutocorrelator& prn_autocorrelator);

    virtual ~PilotedMultiplePrnCorrelator();
    /**
     * Set the PRN source block pointer. It is assumed the samples are one PRN length.
     * \param fftw_source_block source samples
     */
    void set_source_block(wsgc_fftw_complex *fftw_source_block);
    
    /**
     * Execute one correlation operation with the current PRN samples. The Global PRN index is incremented after
     * each call to this method.
     * \param prn_code_index Index of pilot PRN code
     */
    void make_correlation(unsigned int pilot_prn_code_index);
    

protected:
    PilotCorrelationAnalyzer& _pilot_correlation_analyzer; //!< Reference to the pilot correlation analyzer
    std::vector<CorrelationRecord>& _correlation_records; //!< Reference to the message correlation records
    PilotCorrelator& _pilot_correlator; //!< Reference to the pilot correlator
    PilotedMessageCorrelator& _message_correlator; //!< Reference to the message correlator
    PrnAutocorrelator& _prn_autocorrelator; //!< Reference to the PRN autocorrelator
};

#endif // __PILOTED_MULTIPLE_PRN_CORRELATOR_H__
