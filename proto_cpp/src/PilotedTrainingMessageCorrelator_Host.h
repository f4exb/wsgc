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

     PilotedTrainingMessageCorrelator_Host

     Given the frequency and time displacement of correlation peak given by the
     Pilot Correlator it searches correlation for all PRNs of the training sequence symbols
     it shifts accumulates the result to detect a peak corresponding to the PRN at the 
     start of the received sequence. It uses straightforward time correlation.

     This is the host implementation

*/

#ifndef __PILOTED_TRAINING_MESSAGE_CORRELATOR_HOST_H__
#define __PILOTED_TRAINING_MESSAGE_CORRELATOR_HOST_H__

#include "WsgcTypes.h"
#include "ContinuousPhaseCarrier.h"
#include "PilotedTrainingMessageCorrelator.h"
#include <vector>

class LocalCodes_Host;
class PilotCorrelationAnalyzer;

/**
 * \brief Correlator engine to acquire starting point in time of the training sequence using the frequency
 * and time displacement of correlation peak given by the Pilot Correlator - Host implementation
 */
class PilotedTrainingMessageCorrelator_Host : public PilotedTrainingMessageCorrelator
{
public:
    /**
    * Correlator engine for message PRNs
    * \param local_codes PRN signals local copy
    * \param f_sampling Sampling frequency
    * \param f_chip Chip rate (frequency)
    * \param sequence_length Length of training sequence should be less than possible symbol numbers
    */
	PilotedTrainingMessageCorrelator_Host(LocalCodes_Host& local_codes,
			wsgc_float f_sampling,
			wsgc_float f_chip,
			unsigned int sequence_length);
	virtual ~PilotedTrainingMessageCorrelator_Host();

	/**
	 * Do the training sequence correlation over the length of one analysis window.
     * \param pilot_correlation_analyzer Reference to the pilot correlation analyzer
	 */
	virtual void execute(PilotCorrelationAnalyzer& pilot_correlation_analyzer);
    
protected:
    LocalCodes_Host& _local_codes; //!< Reference to the PRN signals local copy.
    ContinuousPhaseCarrier _local_oscillator; //!< Local oscillator for receiving frequency adjustment
    unsigned int _nb_msg_prns; //!< Number of message PRNs in the code set
    unsigned int _fft_N; //!< Size of FFT
    wsgc_complex *_src; //!< Source samples of the current PRN multiplied by local oscillator
    wsgc_complex *_corr_results; //!< Result of correlation of PRNs.
    wsgc_float *_mag_avgsums; //!< Result of PRNs correlation magnitudes shifting averaging sums.
};

#endif /* __PILOTED_TRAINING_MESSAGE_CORRELATOR_HOST_H__ */
