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

     DifferentialModulationMultiplePrnCorrelator

     This flavour of correlator deals with PRNs without use of a pilot sequence
     This is the Host version

*/

#ifndef __UNPILOTED_TRAINING_MESSAGE_CORRELATOR_HOST_H__
#define __UNPILOTED_TRAINING_MESSAGE_CORRELATOR_HOST_H__

#include "UnpilotedTrainingMessageCorrelator.h"
#include "LocalCodesFFT_Host.h"
#include <vector>
#include <cstring>

/**
 * \brief Correlator engine to do correlation of training sequence PRN(s) without use of a pilot sequence. Host version
 */
class UnpilotedTrainingMessageCorrelator_Host : public UnpilotedTrainingMessageCorrelator
{
public:
    /**
    * Correlator engine for message PRNs
    * \param f_sampling Sampling frequency
    * \param f_chip Chip rate (frequency)
    * \param prn_length Length of a PRN sequence in number of chips
    * \param sequence_length Number of PRNs int the training sequence
    * \param averaging_length Number of PRNs used for sliding averaging
    * \param prn_list Reference to the vector of PRN numbers with which to make correlation
    * \param symbol_window_size Number of symbols used for processing. Storage is reserved for symbol_window_size times prn_per_symbol PRN samples
    * \param training_correlation_records Reference to the correlation records
    * \param local_codes_fft_base Reference to the FFT copy of local codes for base modulation
    */
	UnpilotedTrainingMessageCorrelator_Host(
			wsgc_float f_sampling,
			wsgc_float f_chip,
			unsigned int prn_length,
            unsigned int sequence_length,
            unsigned int averaging_length,
			const std::vector<unsigned int>& prn_list,
			std::vector<TrainingCorrelationRecord>& training_correlation_records,
			const LocalCodesFFT_Host& local_codes_fft_base
        );
        
	virtual ~UnpilotedTrainingMessageCorrelator_Host();

	/**
	 * Do the message correlation over the PRN window
     * \param prn_per_symbol Number of PRNs per symbol
	 */
	virtual void execute();

	/**
	 * Append source samples for one PRN length to the buffer
     * \param samples Pointer to source samples
	 */
	virtual void set_samples(wsgc_complex *samples);

protected:
    const LocalCodesFFT_Host& _local_codes_fft_base; //!< Reference to the FFT copy of local codes for base modulation
    wsgc_complex *_samples; //!< Copy of input samples for the length of the window
    wsgc_complex *_src_fft; //!< Samples FFT
    wsgc_complex *_ifft_in_tmp; //!< temporary zone for IFFT input
    wsgc_complex *_ifft_out_tmp; //!< temporary zone for IFFT output
    wsgc_complex *_corr_out_sums; //!< Sum result of correlation of PRNs.
    wsgc_float *_mag_avgsums; //!< Result of PRNs correlation magnitudes shifting averaging sums.
    wsgc_fftw_plan _fft_plan; //!< FFTW plan for forward FFT.
    wsgc_fftw_plan _ifft_plan; //!< FFTW plan for inverse FFT.
    
    static const wsgc_complex _c_zero;

    /**
     * Do a whole correlation process over the PRNs window
     */
    void do_correlation();

    /**
     * Do correlation storing result at the chip sample shift place
     * \param prn_wi PRN index in window
     */
    void do_correlation(unsigned int prn_wi);

    /**
     * Do a sum averaging of correlations over the PRNs window
     * \param prn_wi PRN index in window
     */
    void do_sum_averaging();

};

#endif /* __UNPILOTED_MULTIPLE_PRN_CORRELATOR_HOST_H__ */
