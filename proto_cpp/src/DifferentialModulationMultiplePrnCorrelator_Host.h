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

     This flavour of correlator deals with PRNs encoded with differential modulation
     This is the Host version

*/

#ifndef __DIFFERENTIAL_MODULATION_MULTIPLE_PRN_CORRELATOR_HOST_H__
#define __DIFFERENTIAL_MODULATION_MULTIPLE_PRN_CORRELATOR_HOST_H__

#include "DifferentialModulationMultiplePrnCorrelator.h"
#include "LocalCodesFFT_Host.h"
#include <cstring>

/**
 * \brief Correlator engine to do correlation of PRN(s) encoded with differential modulation. Host version
 */
class DifferentialModulationMultiplePrnCorrelator_Host : public DifferentialModulationMultiplePrnCorrelator
{
public:
    /**
    * Correlator engine for message PRNs
    * \param f_sampling Sampling frequency
    * \param f_chip Chip rate (frequency)
    * \param prn_length Length of a PRN sequence in number of chips
    * \param prn_list Reference to the vector of PRN numbers with which to make correlation
    * \param prn_window_size Number of PRNs used for processing. This is the maximum number of PRNs stored in memory at once
    * \param local_codes_fft_base Reference to the FFT copy of local codes for base modulation
    */
	DifferentialModulationMultiplePrnCorrelator_Host(
			wsgc_float f_sampling,
			wsgc_float f_chip,
			unsigned int prn_length,
			const std::vector<unsigned int>& prn_list,
			unsigned int prn_window_size,
			const LocalCodesFFT_Host& local_codes_fft_base
        );
        
	virtual ~DifferentialModulationMultiplePrnCorrelator_Host();

	/**
	 * Do the message correlation over the PRN window
     * \param prn_per_symbol Number of PRNs per symbol
	 */
	virtual void execute_message();

	/**
	 * Do the training sequence correlation over the PRN window
     * \param prn_per_batch Number of PRNs per sliding average batch
	 */
	virtual void execute_training();
    
	/**
	 * Append source samples for one PRN length to the buffer
     * \param samples Pointer to source samples
	 */
	virtual void set_samples(wsgc_complex *samples);

	/**
	 * Set the initial chip samples for differential processing
     * \param sample Pointer to initial source sample
	 */
	virtual void set_initial_chip(wsgc_complex *chip_samples);

protected:
    const LocalCodesFFT_Host& _local_codes_fft_base; //!< Reference to the FFT copy of local codes for base modulation
    wsgc_complex *_samples; //!< Copy of input samples for the length of the window
    wsgc_complex *_demod; //!< Demodulated samples temporary area
    wsgc_complex *_ifft_in_tmp; //!< temporary zone for IFFT input
    wsgc_complex *_ifft_out_tmp; //!< temporary zone for IFFT output
    wsgc_complex *_corr_out; //!< Correlation results
    wsgc_fftw_plan _ifft_plan; //!< FFTW plan for inverse FFT.
    
    /**
     * Differentially demodulate for the window length by doing a sample to sample multiplication shifted by a chip of samples
     */
    void differentially_demodulate_window();
    
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
     * At end of window processing move last chip to beginning of window
     */
    void chip_carry();
};

#endif /* __DIFFERENTIAL_MODULATION_MULTIPLE_PRN_CORRELATOR_HOST_H__ */
