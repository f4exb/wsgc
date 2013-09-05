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
#ifndef __SINGLE_PRN_CORRELATOR_FREQ_DEP_HOST_H__
#define __SINGLE_PRN_CORRELATOR_FREQ_DEP_HOST_H__

#include "SinglePrnCorrelator_FreqDep.h"
#include "LocalCodesFFT_Host.h"
#include "WsgcTypes.h"
#include <vector>

class GoldCodeGenerator;

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
 * This is a CPU (not GPU) implementation and is intended to be used with the PilotCorrelator and PilotedMultiplePrnCorrelator_Host classes.
 *
 */
class SinglePrnCorrelator_FreqDep_Host : public SinglePrnCorrelator_FreqDep
{
public:
	/**
	* Runs on a source samples batch (one PRN length)
	* \param gc_generator Gold Code generator used to build the codes
	* \param code_modulator Modulator used to build the codes
	* \param f_sampling Sampling frequency
	* \param f_chip Chip rate (frequency)
    * \param _pilot_symbols Reference to the list of pilot symbol PRNs
	* \param nb_f_bins Number of frequency bins explored around IF=0
	* \param prn_per_block Number of PRNs per (symbol) block
	* \param nb_batch_prns Number of PRNs processed in one batch ("PRN batch factor")
	* \param frequency_step_division Frequency step division
	*/
	SinglePrnCorrelator_FreqDep_Host(
			const GoldCodeGenerator& gc_generator,
			CodeModulator& code_modulator,
			wsgc_float f_sampling,
			wsgc_float f_chip,
			std::vector<unsigned int>& pilot_symbols,
			unsigned int nb_f_bins,
			unsigned int prn_per_block=4,
			unsigned int nb_batch_prns=3,
			unsigned int frequency_step_division=1);

	virtual ~SinglePrnCorrelator_FreqDep_Host();

	/**
	 * Runs multiplication and IFFT on a source samples one PRN length.
	 * \param source_block Pointer to first element in FFT of source samples array
	 * \param pilot_prn_index Index of pilot PRN in the local codes matrix
	 * \param PRN position in a two batch (storage depth) cycle
	 */
	void multiply_and_ifft(const wsgc_complex *source_block, unsigned int pilot_prn_index, unsigned int prn_position);

	/**
	 * Do the averaging sum of one half of the IFFT frames
	 * \param first_half true if the first half is to be processed else false
	 */
	virtual void execute_averaging(bool first_half);

	/**
	 * Get the IFFT averaged frames
	 * \return pointer to the beginning of IFFT averaged frames. Time first i.e index by (_nb_f_bins * fi + ti)*storage_depth + bi.
	 */
	const wsgc_complex *get_averaged_ifft_frames() const;


protected:
	LocalCodesFFT_Host _local_codes;  //!< PRN signals local copy. These are the conjugated FFTs in fact.
	wsgc_fftw_plan _ifft_plan;        //!< FFTW plan for inverse FFT.
	wsgc_complex *_ifft_code_in;      //!< Input samples for inverse FFT
	wsgc_complex *_ifft_code_out;     //!< Inverse FFT output samples. This is the result of correlation.
	wsgc_complex *_ifft_code_in_tmp;  //!< Input samples for inverse FFT temporary storage
	wsgc_complex *_ifft_code_out_tmp; //!< Output samples for inverse FFT temporary storage
};

#endif // __SINGLE_PRN_CORRELATOR_FREQ_DEP_HOST_H__
