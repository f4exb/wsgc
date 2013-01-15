/*
 * PilotedTrainingMultiplePrnCorrelator.h
     Copyright 2012-2013 Edouard Griffiths <f4exb at free dot fr>

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

     PilotedTrainingMultiplePrnCorrelator

     Processes the correlation of training PRNs using the training identifying pilot PRN
     It works with a pilot correlator and a training message correlator referenced from outside.
     - The pilot correlator looks for a correlation on pilot PRN(s) only and produces the information
       on the best possible frequency shift from zero IF and time delta of the PRN origin in the receiver time reference
     - The training message correlator uses this information to make a correlation at the time delay and center frequency
       given by the pilot for all possible PRNs. It accumulates the correlation peaks of all PRNs with appropriate 
       shift so that the sequence of training PRNs would make the peak accumulation on the same PRN bin. Normally one 
       would use an incremental sequence of PRNs (by index order) and shift results one place to the right each time.

 */

#ifndef __PILOTED_TRAINING_MULTIPLE_PRN_CORRELATOR_H__
#define __PILOTED_TRAINING_MULTIPLE_PRN_CORRELATOR_H__

#include "WsgcTypes.h"
#include <vector>
#include <sstream>
#include <time.h>

class PilotCorrelationAnalyzer;
class PilotCorrelator;
class PilotedTrainingMessageCorrelator;

/**
 * \brief Correlator engine to get the message epoch estimation of PRNs using a training sequence. It uses a pilot correlator and a training message correlator
 *
 *  It works with a pilot correlator and a training message correlator referenced from outside.
 *  - The pilot correlator looks for a correlation on pilot PRN(s) only and produces the information
 *    on the best possible frequency shift from zero IF and time delta of the PRN origin in the receiver time reference
 *  - The training message correlator uses this information to make a correlation at the time delay and center frequency
 *    given by the pilot for all possible PRNs. It accumulates the correlation peaks of all PRNs with appropriate 
 *    shift so that the sequence of training PRNs would make the peak accumulation on the same PRN bin. Normally one 
 *    would use an incremental sequence of PRNs (by index order) and shift results one place to the right each time.
 */
class PilotedTrainingMultiplePrnCorrelator
{
public:
    
    /**
    * Training correlator engine constructor
    * \param pilot_correlation_analyzer Reference to the pilot correlation analyzer
    * \param pilot_correlator Reference to the pilot correlator
    * \param training_message_correlator Reference to the training message correlator
    */
    PilotedTrainingMultiplePrnCorrelator(
            PilotCorrelationAnalyzer& pilot_correlation_analyzer,
            PilotCorrelator& pilot_correlator, 
            PilotedTrainingMessageCorrelator& message_correlator);

    virtual ~PilotedTrainingMultiplePrnCorrelator();
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
    PilotCorrelator& _pilot_correlator; //!< Reference to the pilot correlator
    PilotedTrainingMessageCorrelator& _training_message_correlator; //!< Reference to the training message correlator
};

#endif // __PILOTED_TRAINING_MULTIPLE_PRN_CORRELATOR_H__
