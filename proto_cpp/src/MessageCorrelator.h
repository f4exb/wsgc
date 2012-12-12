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

#ifndef __MESSAGE_CORRELATOR_H__
#define __MESSAGE_CORRELATOR_H__

#include "WsgcTypes.h"
#include "CorrelationRecord.h"
#include "ContinuousPhaseCarrier.h"
#include "MessageCorrelationMatrices.h"
#include <vector>

class LocalCodes_Host;
class PilotCorrelationAnalyzer;

/**
 * \brief Correlator engine to acquire and track message PRN(s) using the frequency
 * and time displacement of correlation peak given by the Pilot Correlator
 */
class MessageCorrelator
{
public:
    /**
    * Correlator engine for message PRNs
    * \param local_codes PRN signals local copy
    * \param f_sampling Sampling frequency
    * \param f_chip Chip rate (frequency)
    * \param prn_per_symbol Number of PRNs per symbol or averaging block
    */
	MessageCorrelator(LocalCodes_Host& local_codes, wsgc_float f_sampling, wsgc_float f_chip, unsigned int prn_per_symbol);
	~MessageCorrelator();

	/**
	 * Do the message correlation over the length of one analysis window.
     * \param pilot_correlation_analyzer Reference to the pilot correlation analyzer
	 */
	void execute(PilotCorrelationAnalyzer& pilot_correlation_analyzer);
    
protected:
    LocalCodes_Host& _local_codes; //!< Reference to the PRN signals local copy.
    ContinuousPhaseCarrier _local_oscillator; //!< Local oscillator for receiving frequency adjustment
    wsgc_float _f_sampling; //!< Sampling frequency
    wsgc_float _f_chip; //!< Chip rate
    unsigned int _nb_msg_prns; //!< Number of message PRNs to explore
    unsigned int _prn_per_symbol; //!< Number of PRNs per symbol that is also the number of PRNs per averaging block
    unsigned int _fft_N; //!< Size of FFT
    wsgc_float _delta_f; //!< Retain receiving frequency
    wsgc_complex *_src; //!< Source samples of the current PRN multiplied by local oscillator
    wsgc_complex *_corr_results; //!< Result of correlation of PRNs.
    wsgc_complex *_corr_avgsums; //!< Result of PRNs correlation averages. Sparse matrix PRNi x [(delta t, max)]times PRN per symbol
    wsgc_complex *_noise_corr_results; //!< Result of correlation of noise PRN.
};

#endif /* __MESSAGE_CORRELATOR_H__ */
