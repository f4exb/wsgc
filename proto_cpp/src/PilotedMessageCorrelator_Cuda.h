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

     MessageCorrelator_Host

     Given the frequency and time displacement of correlation peak given by the
     Pilot Correlator it searches correlation for all PRNs of the message symbols
     to select the one that was sent. It uses straightforward time correlation.

     This is the CUDA implementation

*/

#ifndef __MESSAGE_CORRELATOR_CUDA_H__
#define __MESSAGE_CORRELATOR_CUDA_H__

#include "WsgcTypes.h"
#include "CorrelationRecord.h"
#include "ContinuousPhaseCarrier.h"
#include "PilotedMessageCorrelator.h"

#include <thrust/device_vector.h>
#include <cuComplex.h>

#include <vector>

class LocalCodes_Cuda;
class PilotCorrelationAnalyzer;
// TODO: make the CUDA and Host classes create their own flavour of local codes object. Pass Gold Code generator and modulator references to the super class.

/**
 * \brief Correlator engine to acquire and track message PRN(s) using the frequency
 * and time displacement of correlation peak given by the Pilot Correlator - Host implementation
 */
class PilotedMessageCorrelator_Cuda : public PilotedMessageCorrelator
{
public:
    /**
    * Correlator engine for message PRNs
    * \param local_codes PRN signals local copy
    * \param f_sampling Sampling frequency
    * \param f_chip Chip rate (frequency)
    * \param prn_per_symbol Number of PRNs per symbol or averaging block
    */
	PilotedMessageCorrelator_Cuda(LocalCodes_Cuda& local_codes, wsgc_float f_sampling, wsgc_float f_chip, unsigned int prn_per_symbol);
	virtual ~PilotedMessageCorrelator_Cuda();

	/**
	 * Do the message correlation over the length of one analysis window.
     * \param pilot_correlation_analyzer Reference to the pilot correlation analyzer
	 */
	virtual void execute(PilotCorrelationAnalyzer& pilot_correlation_analyzer);

protected:
    LocalCodes_Cuda& _local_codes; //!< Reference to the PRN signals local copy.
    ContinuousPhaseCarrier _local_oscillator; //!< Local oscillator for receiving frequency adjustment
    unsigned int _nb_msg_prns; //!< Number of message PRNs to explore
    unsigned int _fft_N; //!< Size of FFT
    wsgc_complex *_src; //!< Source samples of the current PRN multiplied by local oscillator
    thrust::device_vector<cuComplex> _d_corr_in; //!< Frequency mixed source to correlate
    thrust::device_vector<cuComplex> _d_mul_out; //!< Result of fixed delay multiplications, one FFT size block per PRN, stride for averaging (average first strategy)
    thrust::device_vector<cuComplex> _d_corr_out; //!< Correlation result (reduce by key values), one sample per PRN
    thrust::device_vector<cuComplex> _d_corr_out_avg; //!< intermediate results for averaging
    thrust::device_vector<int> _d_keys; //!< Keys for the reduce by key (unused), one key per PRN

};

#endif /* __MESSAGE_CORRELATOR_CUDA_H__ */
