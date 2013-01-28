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

*/

#ifndef __DIFFERENTIAL_MODULATION_MULTIPLE_PRN_CORRELATOR_H__
#define __DIFFERENTIAL_MODULATION_MULTIPLE_PRN_CORRELATOR_H__

#include "WsgcTypes.h"
#include "CorrelationRecord.h"
#include <vector>

/**
 * \brief Correlator engine to do correlation of PRN(s) encoded with differential modulation
 */
class DifferentialModulationMultiplePrnCorrelator
{
public:
    /**
    * Correlator engine for message PRNs
    * \param f_sampling Sampling frequency
    * \param f_chip Chip rate (frequency)
    * \param prn_length Length of a PRN sequence in number of chips
    * \param prn_list Reference to the vector of PRN numbers with which to make correlation
    * \param prn_window_size Number of PRNs used for processing. This is the maximum number of PRNs stored in memory at once
    */
	DifferentialModulationMultiplePrnCorrelator(
			wsgc_float f_sampling,
			wsgc_float f_chip,
			unsigned int prn_length,
			const std::vector<unsigned int>& prn_list,
			unsigned int prn_window_size);
        
	virtual ~DifferentialModulationMultiplePrnCorrelator();

	/**
	 * Do the message correlation over the PRN window
     * \param prn_per_symbol Number of PRNs per symbol
	 */
	virtual void execute_message() = 0;

	/**
	 * Do the training sequence correlation over the PRN window
     * \param prn_per_batch Number of PRNs per sliding average batch
	 */
	virtual void execute_training() = 0;

	/**
	 * Append source samples for one PRN length to the buffer
     * \param samples Pointer to source samples
	 */
	virtual void set_samples(wsgc_complex *samples) = 0;

	/**
	 * Set the initial chip samples for differential processing
     * \param chip_samples Pointer to initial source samples
	 */
	virtual void set_initial_chip(wsgc_complex *chip_samples) = 0;
    
protected:
    wsgc_float _f_sampling; //!< Sampling frequency
    wsgc_float _f_chip; //!< Chip rate
    unsigned int _prn_length; //!< Length of a PRN sequence
    unsigned int _fft_N; //!< Length of FFT that is also the length of a PRN sequence in number of samples
    wsgc_float _fractional_chip_per_sample; //!< Length of one sample in chip unit (should be < 1/4)
    wsgc_float _fractional_samples_per_chip; //!< Fractional number of samples per chip (should be > 4)
    unsigned int _int_samples_per_chip; //!< Integral number of samples per chip (floor of above)
    const std::vector<unsigned int>& _prn_list; //!< Reference to the vector of PRN numbers with which to make correlation
    unsigned int _prn_window_size; //!<  Number of PRNs used for processing one execution (single call to "execute"). This is the maximum number of PRNs stored in memory
    unsigned int _samples_length; //!< Number of stored source samples
};

#endif /* __DIFFERENTIAL_MODULATION_MULTIPLE_PRN_CORRELATOR_H__ */
