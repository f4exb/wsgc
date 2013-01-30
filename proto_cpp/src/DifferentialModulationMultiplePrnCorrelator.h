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
#include <vector>

class TrainingCorrelationRecord;
class CorrelationRecord;

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
    * \param prn_per_symbol Number of PRNs per symbol
    * \param prn_list Reference to the vector of PRN numbers with which to make correlation
    * \param symbol_window_size Number of symbols used for processing. Storage is reserved for symbol_window_size times prn_per_symbol PRN samples
    * \param correlation_records Reference to the correlation records
    * \param training_correlation_records Reference to the training correlation records
    */
	DifferentialModulationMultiplePrnCorrelator(
			wsgc_float f_sampling,
			wsgc_float f_chip,
			unsigned int prn_length,
            unsigned int prn_per_symbol,
			const std::vector<unsigned int>& prn_list,
			unsigned int symbol_window_size,
			std::vector<CorrelationRecord>& correlation_records,
			std::vector<TrainingCorrelationRecord>& training_correlation_records
			);
        
	virtual ~DifferentialModulationMultiplePrnCorrelator();

	/**
	 * Do the message correlation over the symbols window.
	 */
	virtual void execute_message() = 0;

	/**
	 * Do the training sequence correlation over the symbols window.
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

	/**
	 * Tells if the window is ready for processing
	 */
	bool is_window_ready()
	{
		return _prns_length == _symbol_window_size * _prn_per_symbol;
	}
    
protected:
    wsgc_float _f_sampling; //!< Sampling frequency
    wsgc_float _f_chip; //!< Chip rate
    unsigned int _prn_length; //!< Length of a PRN sequence
    unsigned int _fft_N; //!< Length of FFT that is also the length of a PRN sequence in number of samples
    unsigned int _prn_per_symbol; //!< Number of PRNs per message symbol
    unsigned int _global_prn_index; //!< Index of PRN in all received samples
    wsgc_float _fractional_chip_per_sample; //!< Length of one sample in chip unit (should be < 1/4)
    wsgc_float _fractional_samples_per_chip; //!< Fractional number of samples per chip (should be > 4)
    unsigned int _int_samples_per_chip; //!< Integral number of samples per chip (floor of above)
    const std::vector<unsigned int>& _prn_list; //!< Reference to the vector of PRN numbers with which to make correlation
    unsigned int _symbol_window_size; //!<  Window of symbols in use for processing.
    unsigned int _samples_length; //!< Number of stored source samples
    unsigned int _prns_length; //!< Number of stored PRNs
    std::vector<wsgc_float> _max_sy_mags; //!< Maximum magnitudes for each symbol
    std::vector<wsgc_float> _max_sy_mags_prni; //!< Maximum magnitudes sums (at each PRN) for each symbol
    std::vector<unsigned int> _max_sy_iffti; //!< Maximum ifft indexes at maximum magnitudes for each symbol
    std::vector<unsigned int> _max_sy_prni; //!< Maximum ifft indexes at maximum magnitudes for each symbol
	std::vector<CorrelationRecord>& _correlation_records; //!< Reference to the correlation records
	std::vector<TrainingCorrelationRecord>& _training_correlation_records; //!< Reference to the training correlation records


    void init_results();
};

#endif /* __DIFFERENTIAL_MODULATION_MULTIPLE_PRN_CORRELATOR_H__ */
