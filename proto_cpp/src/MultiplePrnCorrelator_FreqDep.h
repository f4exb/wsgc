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
      
     MultiplePrnCorrelator_FreqDep
      
     Does a IFFT batch correlation and averaging given a list of PRNs to look for over a range of frequencies.
     
*/
#ifndef __MULTIPLE_PRN_CORRELATOR_FREQ_DEP_H__
#define __MULTIPLE_PRN_CORRELATOR_FREQ_DEP_H__

#include "WsgcTypes.h"
#include <vector>

class GoldCodeGenerator;
class LocalCodesFFT;

/**
 * \brief Correlator engine to look for a list of PRN(s) in a range of shift frequencies. The list of PRNs is given by the Local Codes object.
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
class MultiplePrnCorrelator_FreqDep
{
public:
	/**
	* Runs on a source samples batch (one PRN length)
	* \param fft_N FFT size
	* \param nb_f_bins Number of frequency bins explored around IF=0
	* \param prn_per_symbol Number of PRNs per (symbol) block
	* \param nb_batch_prns Number of PRNs processed in one batch ("PRN batch factor")
	* \param frequency_step_division Frequency step division
	*/
	MultiplePrnCorrelator_FreqDep(
			unsigned int fft_N,
			wsgc_float f_sampling,
			unsigned int nb_message_symbols,
			unsigned int nb_f_bins,
			unsigned int prn_per_symbol=4,
			unsigned int nb_batch_prns=3,
			unsigned int frequency_step_division=1);

	virtual ~MultiplePrnCorrelator_FreqDep();
    
	/**
	 * Runs multiplication and IFFT on a source samples one PRN length.
	 * \param source_block Pointer to first element in FFT of source samples array
	 * \param PRN position in a two batch (storage depth) cycle
	 */
    virtual void multiply_and_ifft(const wsgc_complex *source_block, unsigned int prn_position) = 0;
    
	/**
	 * Do the averaging sum of one half of the IFFT frames
	 * \param first_half true if the first half is to be processed else false
	 */
	virtual void execute_averaging(bool first_half) = 0;

	/**
	 * Get the vector of maximum magnitudes for this batch
     * \param prni Index of PRNs in the list of PRNs to look for
	 * \return reference to the maximum magnitudes for this batch vector
	 */
	const std::vector<wsgc_float>& get_batch_magnitudes_max(unsigned int prni) const
	{
		return _batch_max_magnitudes[prni];
	}

	/**
	 * Get the vector of composite indexes of magnitude maxima (_nb_f_bins * (_frequency_step_division * fi + fsi) + ti) for this batch
     * \param prni Index of PRNs in the list of PRNs to look for
	 * \return reference to the composite indexes of magnitude maxima vector
	 */
	const std::vector<unsigned int>& get_batch_composite_indexes_max(unsigned int prni) const
	{
		return _batch_max_composite_indexes[prni];
	}

	/**
	 * Get the vector of complex values at magnitude maxima for this batch
     * \param prni Index of PRNs in the list of PRNs to look for
	 * \return reference to the complex values at magnitude maxima vector
	 */
	const std::vector<wsgc_complex>& get_batch_complex_values_max(unsigned int prni) const
	{
		return _batch_complex_values_max[prni];
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
    
    /**
     * Get the number of PRNs processed in one batch
     */
    unsigned int get_nb_batch_prns() const
    {
        return _nb_batch_prns;
    }
    
    /**
     * Get the storage depth factor: 2*(nb_batch_prns)
     */
    unsigned int get_storage_depth() const
    {
        return _storage_depth;
    }
    
    /**
     * Get the number of PRNs (codes) explored
     */
    unsigned int get_nb_codes() const
    {
        return _nb_codes;
    }
    
    /**
     * Get the current batch index (0 first)
     */
    unsigned int get_batch_index() const
    {
        return _batch_index;
    }

    /**
     * Get the number of PRNs per symbol
     */
    unsigned int get_prn_per_symbol() const
    {
    	return _prn_per_symbol;
    }

    /**
     * Get the FFT size
     */
    unsigned int get_fft_N() const
    {
    	return _fft_N;
    }

    /**
     * Get the frequency step division
     */
    unsigned int get_freq_step_division() const
    {
    	return _freq_step_division;
    }

    /**
     * Get the sampling frequency
     */
    wsgc_float get_f_sampling() const
    {
    	return _f_sampling;
    }

    /**
     * Get the number of frequency bins explored around IF=0
     */
    unsigned int get_nb_f_bins() const
    {
    	return _nb_f_bins;
    }

    /**
     * Get the number of possible message symbols
     */
    unsigned int get_nb_message_symbols() const
    {
        return _nb_message_symbols;
    }

protected:
	unsigned int _fft_N;               //!< FFT size
	wsgc_float    _f_sampling;         //!< Sampling frequency
	unsigned int _nb_message_symbols;  //!< Number of possible message symbols
	unsigned int _nb_f_bins;           //!< Number of frequency bins explored around IF=0
	unsigned int _freq_step_division;  //!< Further division of the frequency step. Corresponds to the number of FFT inputs
	unsigned int _prn_per_symbol;      //!< Number of PRNs per symbol
	unsigned int _nb_batch_prns;       //!< Number of PRNs processed in one batch
    unsigned int _nb_codes;            //!< Number of PRNs to look for. These are the number of different codes stored in the Local Codes object.

	unsigned int _storage_depth;       //!< Storage depth factor: 2*(nb_batch_prns)
	unsigned int _batch_index;         //!< Batch index. Incremented after each batch processing. Ever increasing.

	std::vector<std::vector<wsgc_float> > _batch_max_magnitudes;          //!< Max magnitudes for the current batch (one per code)
	std::vector<std::vector<unsigned int> > _batch_max_composite_indexes; //!< Corresponding composite indexes (_nb_f_bins * fi + ti) (one per code)
	std::vector<std::vector<wsgc_complex> > _batch_complex_values_max;    //!< Max complex values for the current batch (one per code)

	static const wsgc_complex cone; //!< One complex value

};

#endif // __MULTIPLE_PRN_CORRELATOR_FREQ_DEP_H__
