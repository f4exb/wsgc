/*
     Copyright 2012 Edouard Griffiths <f4exb at free dot fr>
 
     This file is part of WSGC. A Weak Signal transmission mode using Gold Codes.
 
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
      
     SinglePrnCorrelator_FreqDep
      
     Does a IFFT batch correlation and averaging given one single PRN to look for over a range of frequencies.
     
*/
#ifndef __SINGLE_PRN_CORRELATOR_FREQ_DEP_H__
#define __SINGLE_PRN_CORRELATOR_FREQ_DEP_H__

#include "WsgcTypes.h"
#include <vector>

class GoldCodeGenerator;
class LocalCodesFFT;

/**
 * \brief Correlator engine to look for one PRN(s) in a range of shift frequencies
 *
 * Original one PRN length samples have been multiplied by a range of sub-frequencies in a previous step
 * Takes the succession of sub-frequency mixed PRN samples and processes it in order to:
 * - process correlation using multiplication in the frequency domain and IFFT method
 * - do the complex averaging 
 *
 * Does correlation in the frequency domain using multiplication by the local code (times the number of frequency shift steps) followed by IFFT.
 * Frequency shift is performed by rotating the local code. This gives discrete frequencies which step corresponds to the FFT bin size that is
 * also the recurrence frequency of the PRN code. This leaves gaps where correlation cannot be made as it has to be at most 10% of the FFT bin size
 * away. Therefore at each step the bunch of pre-mixed PRN copies is processed (and not only just one PRN at zero IF). The pre-mixing is done at
 * an earlier stage before this correlator is invoked.
 * In addition it takes a block of inputs to be processed in one batch (see prn_per_block parameter of constructor).
 *
 * This is the common abstract class for CPU and  GPU implementations
 *
 */
class SinglePrnCorrelator_FreqDep
{
public:
	/**
	* Runs on a source samples batch (one PRN length)
	* \param fft_N FFT size
	* \param nb_f_bins Number of frequency bins explored around IF=0
	* \param prn_per_block Number of PRNs per (symbol) block
	* \param nb_batch_prns Number of PRNs processed in one batch ("PRN batch factor")
	* \param frequency_step_division Frequency step division
	*/
	SinglePrnCorrelator_FreqDep(
			unsigned int fft_N,
			unsigned int nb_f_bins,
			unsigned int prn_per_block=4,
			unsigned int nb_batch_prns=3,
			unsigned int frequency_step_division=1);

	virtual ~SinglePrnCorrelator_FreqDep();

	/**
	 * Do the averaging sum of one half of the IFFT frames
	 * \param first_half true if the first half is to be processed else false
	 */
	virtual void execute_averaging(bool first_half) = 0;

	/**
	 * Get the vector of maximum magnitudes for this batch
	 * \return reference to the maximum magnitudes for this batch vector
	 */
	const std::vector<wsgc_float>& get_batch_magnitudes_max() const
	{
		return _batch_max_magnitudes;
	}

	/**
	 * Get the vector of composite indexes of magnitude maxima (_nb_f_bins * (_frequency_step_division * fi + fsi) + ti) for this batch
	 * \return reference to the composite indexes of magnitude maxima vector
	 */
	const std::vector<unsigned int>& get_batch_composite_indexes_max() const
	{
		return _batch_max_composite_indexes;
	}

	/**
	 * Get the vector of complex values at magnitude maxima for this batch
	 * \return reference to the complex values at magnitude maxima vector
	 */
	const std::vector<wsgc_complex>& get_batch_complex_values_max() const
	{
		return _batch_complex_values_max;
	}

	/**
	 * Calculate individual indexes from a composite index value
	 * \param composite_index Composite index
	 * \param f_index Reference to the frequency bin that will be updated
	 * \param fs_index Reference to the frequency sub-step that will be updated
	 * \param t_index Reference to the time shift index that will be updated
	 */
	void calculate_indexes(unsigned int composite_index, unsigned int& f_index, unsigned int& fs_index, unsigned int& t_index)
	{
		t_index = composite_index % _fft_N;
		fs_index = (composite_index / _fft_N) % _freq_step_division;
		f_index = (composite_index / _fft_N) / _freq_step_division;
	}

protected:
	unsigned int _fft_N;                       //!< FFT size
	unsigned int _nb_f_bins;                   //!< Number of frequency bins explored around IF=0
	unsigned int _freq_step_division;          //!< Further division of the frequency step. Corresponds to the number of FFT inputs
	unsigned int _prn_per_block;               //!< Number of PRNs per averaging block
	unsigned int _nb_batch_prns;               //!< Number of PRNs processed in one batch

	unsigned int _storage_depth;       //!< Storage depth factor: 2*(nb_batch_prns)
	unsigned int _batch_index;         //!< Batch index. Incremented after each batch processing. Ever increasing.

	std::vector<wsgc_float> _batch_max_magnitudes;          //!< Max magnitudes for the current batch
	std::vector<unsigned int> _batch_max_composite_indexes; //!< Corresponding composite indexes (_nb_f_bins * fi + ti)
	std::vector<wsgc_complex> _batch_complex_values_max;    //!< Max complex values for the current batch;

	static const wsgc_complex cone; //!< One complex value

};

#endif // __SINGLE_PRN_CORRELATOR_FREQ_DEP_H__
