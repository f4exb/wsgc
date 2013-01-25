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

     PilotedTrainingMessageCorrelator_Cuda

     Given the frequency and time displacement of correlation peak given by the
     Pilot Correlator it searches correlation for all PRNs of the training sequence symbols
     it shifts accumulates the result to detect a peak corresponding to the PRN at the 
     start of the received sequence. It uses straightforward time correlation.

     This is the CUDA implementation

*/

#ifndef __PILOTED_TRAINING_MESSAGE_CORRELATOR_CUDA_H__
#define __PILOTED_TRAINING_MESSAGE_CORRELATOR_CUDA_H__

#include "WsgcTypes.h"
#include "ContinuousPhaseCarrier.h"
#include "PilotedTrainingMessageCorrelator.h"

#include <thrust/device_vector.h>
#include <cuComplex.h>

#include <vector>

class LocalCodes_Cuda;
class PilotCorrelationAnalyzer;

/**
 * \brief Correlator engine to acquire starting point in time of the training sequence using the frequency
 * and time displacement of correlation peak given by the Pilot Correlator - Host implementation
 */
class PilotedTrainingMessageCorrelator_Cuda : public PilotedTrainingMessageCorrelator
{
public:
    /**
    * Correlator engine for message PRNs
    * \param local_codes PRN signals local copy
    * \param f_sampling Sampling frequency
    * \param f_chip Chip rate (frequency)
    * \param analysis_window_size Analysis window size in number of PRNs
    * \param sequence_length Length of training sequence should be less than possible symbol numbers
    */
	PilotedTrainingMessageCorrelator_Cuda(LocalCodes_Cuda& local_codes, 
            wsgc_float f_sampling, 
            wsgc_float f_chip, 
            unsigned int analysis_window_size,
            unsigned int sequence_length);
            
	virtual ~PilotedTrainingMessageCorrelator_Cuda();

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
    unsigned int _analysis_window_size; //!< Analysis window size in number of PRNs
    wsgc_complex *_src; //!< Source samples of the current PRN multiplied by local oscillator
    thrust::device_vector<cuComplex> _d_corr_in; //!< Frequency mixed source to correlate
    thrust::device_vector<cuComplex> _d_mul_out; //!< Result of fixed delay multiplications, one FFT size block per PRN, stride for averaging (average first strategy)
    thrust::device_vector<cuComplex> _d_corr_out; //!< Correlation result (reduce by key values), one sample per PRN
    thrust::device_vector<float> _d_corr_mag; //!< Correlation result magnitudes. Magnitude (squared norm) of the above
    thrust::device_vector<float> _d_corr_mag_avgsum; //!< averaging sum results
    thrust::device_vector<int> _d_keys; //!< Keys for the reduce by key (unused), one key per PRN
    thrust::host_vector<int> _h_keys; //!< Keys for the reduce by key (unused), one key per PRN
    thrust::host_vector<float> _h_corr_mag_avgsum; //!< Correlation result magnitudes (host copy)
    static const cuComplex _c_zero; //!< Complex zero for initializations
    wsgc_float _max_avg; //!< Maximum shifted average sum of magnitudes
    unsigned int _max_avg_index ; //!< Index of the maximum
};

#endif /* __PILOTED_TRAINING_MESSAGE_CORRELATOR_CUDA_H__ */
