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

     MessageCorrelator

     Given the frequency and time displacement of correlation peak given by the
     Pilot Correlator it searches correlation for all PRNs of the message symbols
     to select the one that was sent. It uses straightforward time correlation.

*/

#ifndef __PILOTED_MESSAGE_CORRELATOR_H__
#define __PILOTED_MESSAGE_CORRELATOR_H__

#include "WsgcTypes.h"
#include "CorrelationRecord.h"
#include "ContinuousPhaseCarrier.h"
#include <vector>

class PilotCorrelationAnalyzer;

/**
 * \brief Correlator engine to acquire and track message PRN(s) using the frequency
 * and time displacement of correlation peak given by the Pilot Correlator
 */
class PilotedMessageCorrelator
{
public:
    /**
    * Correlator engine for message PRNs
    * \param f_sampling Sampling frequency
    * \param f_chip Chip rate (frequency)
    * \param prn_per_symbol Number of PRNs per symbol or averaging block
    */
	PilotedMessageCorrelator(wsgc_float f_sampling, wsgc_float f_chip, unsigned int prn_per_symbol);
	virtual ~PilotedMessageCorrelator();

	/**
	 * Do the message correlation over the length of one analysis window.
     * \param pilot_correlation_analyzer Reference to the pilot correlation analyzer
	 */
	virtual void execute(PilotCorrelationAnalyzer& pilot_correlation_analyzer) = 0;
    
protected:
    wsgc_float _f_sampling; //!< Sampling frequency
    wsgc_float _f_chip; //!< Chip rate
    unsigned int _prn_per_symbol; //!< Number of PRNs per symbol that is also the number of PRNs per averaging block
    wsgc_float _delta_f; //!< Retain receiving frequency
};

#endif /* __PILOTED_MESSAGE_CORRELATOR_H__ */
