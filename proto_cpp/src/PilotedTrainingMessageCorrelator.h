/*
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

     PilotedTrainingMessageCorrelator

     Given the frequency and time displacement of correlation peak given by the
     Pilot Correlator it searches correlation for all PRNs of the training message symbols
     and accumulates results in order to find a peak of magnitude among PRNs.
     The training sequence is the sucession of message PRNs in their PRN index order for the 
     given code set. Accumulation is shifted one PRN place to the right each time (that is one 
     PRN earlier in the index sequence) so that peaks would accumulate always in the same place.
     The PRN index eventually selected gives the relative position in the whole sequence when 
     the process was started. This reveals the point in time when the first PRN in the sequence 
     was sent hence the epoch of the message.
     It uses straightforward time correlation.

*/

#ifndef __PILOTED_TRAINING_MESSAGE_CORRELATOR_H__
#define __PILOTED_TRAINING_MESSAGE_CORRELATOR_H__

#include "WsgcTypes.h"
#include "ContinuousPhaseCarrier.h"
#include <vector>

class PilotCorrelationAnalyzer;

/**
 * \brief Correlator engine to estimate the time epoch of the message
 * Given the frequency and time displacement of correlation peak given by the
 * Pilot Correlator it searches correlation for all PRNs of the training message symbols
 * and accumulates results in order to find a peak of magnitude among PRNs.
 * The training sequence is the sucession of message PRNs in their PRN index order for the 
 * given code set. Accumulation is shifted one PRN place to the right each time (that is one 
 * PRN earlier in the index sequence) so that peaks would accumulate always in the same place.
 * The PRN index eventually selected gives the relative position in the whole sequence when 
 * the process was started. This reveals the point in time when the first PRN in the sequence 
 * was sent hence the epoch of the message.
 * It uses straightforward time correlation.
 */
class PilotedTrainingMessageCorrelator
{
public:
    /**
    * Correlator engine for message PRNs
    * \param f_sampling Sampling frequency
    * \param f_chip Chip rate (frequency)
    * \param sequence_length Length of training sequence should be less than possible symbol numbers
    */
	PilotedTrainingMessageCorrelator(wsgc_float f_sampling, wsgc_float f_chip, unsigned int sequence_length);
	virtual ~PilotedTrainingMessageCorrelator();

	/**
	 * Do the message correlation over the length of one analysis window.
     * \param pilot_correlation_analyzer Reference to the pilot correlation analyzer
	 */
	virtual void execute(PilotCorrelationAnalyzer& pilot_correlation_analyzer) = 0;

    wsgc_float _maxtoavg_max; //!< Maximum ratio between maximum magnitude and the averaging sum of all magnitudes
    wsgc_float _mag_max; //!< Maximum magnitude
    unsigned int _prni_max; //!< PRN index at maximum magnitude
    unsigned int _prnai_max; //!< Analysis window index of best max to average found

protected:
    wsgc_float _f_sampling; //!< Sampling frequency
    wsgc_float _f_chip; //!< Chip rate
    unsigned int _sequence_length; //!< Length of sequence should be less than possible symbol numbers
    wsgc_float _delta_f; //!< Retain receiving frequency
};

#endif /* __PILOTED_TRAINING_MESSAGE_CORRELATOR_H__ */
