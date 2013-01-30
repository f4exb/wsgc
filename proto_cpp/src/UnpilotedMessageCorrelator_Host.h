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

     UnpilotedMessageCorrelator

	 Does correlation on all possible message PRNs (normally also including noise PRN)
	 to find a correlation peak representing the time delay of the start of sequence.

	 Uses frequency domain correlation to check all possible delays at once which is
	 much more computationnaly efficient

	 This is the class for Host implementation

*/

#ifndef __UNPILOTED_MESSAGE_CORRELATOR_HOST_H__
#define __UNPILOTED_MESSAGE_CORRELATOR_HOST_H__

#include "WsgcTypes.h"
#include "UnpilotedMessageCorrelator.h"
#include "LocalCodesFFT_Host.h"
#include <vector>

class CorrelationRecord;
class CodeModulator;
class GoldCodeGenerator;

/**
 * \brief Correlator engine on all possible message PRNs (normally also including noise PRN)
 * to find a correlation peak representing the time delay of the start of sequence.
 *
 * Uses frequency domain correlation to check all possible delays at once which is
 * much more computationnaly efficient
 *
 * This is the class for Host implementation
 */
class UnpilotedMessageCorrelator_Host : public UnpilotedMessageCorrelator
{
public:
    /**
    * Correlator engine for message PRNs
    * \param f_sampling Sampling frequency
    * \param f_chip Chip rate (frequency)
    * \param simulate_symbol_synchronization True if message is externally synchronized therefore the start PRN in symbol sequence is known
    * \param prn_per_symbol Number of PRNs per symbol or averaging block
    * \param nb_batch_prns Number of PRNs per process batch. If synchro is active this is multiplied by two times the Number of PRNs per symbol
    * \param message_symbols List of message symbols to correlate
    * \param code_modulator Code modulator to be used for local signal generation
    * \param gc_generator Gold Code generator to be used for local signal generation
    */
	UnpilotedMessageCorrelator_Host(
			wsgc_float f_sampling,
			wsgc_float f_chip,
			bool simulate_symbol_synchronization,
			unsigned int prn_per_symbol,
			unsigned int nb_batch_prns,
			std::vector<unsigned int>& message_symbols,
			CodeModulator& code_modulator,
			GoldCodeGenerator& gc_generator
			);

	virtual ~UnpilotedMessageCorrelator_Host();

	/**
	 * Assign one PRN length of samples to be processed
	 * \param source_block Samples the length of one PRN to be processed
	 */
	virtual void set_source_block(wsgc_complex *source_block);

	/**
	 * Do the message correlation over the length of one analysis window.
	 * \param message_correlation_record Reference to the message correlation records
     * \return True if a new averaging took plce hence new results are available
	 */
	virtual bool execute(std::vector<CorrelationRecord>& message_correlation_records);
    
protected:
	LocalCodesFFT_Host _local_codes;  //!< Local copy of the FFT conjugate of codes to be correlated
	unsigned int _fft_N;              //!< Size of FFT/IFFT. This is also the number of samples per PRN
	unsigned int _nb_msg_prns;        //!< Number of message PRNs including noise PRN
    wsgc_fftw_plan _fft_sample_plan;  //!< FFTW plan for source FFT
    wsgc_complex *_fft_sample_in;     //!< Input samples
    wsgc_complex *_fft_sample_out;    //!< FFT of input samples
	wsgc_fftw_plan _ifft_plan;        //!< FFTW plan for inverse FFT.
	wsgc_complex *_ifft_code_out;     //!< Inverse FFT output samples. This is the result of correlation.
	wsgc_float   *_ifft_out_mags;     //!< Inverse FFT output samples magnitudes.
	wsgc_complex *_ifft_code_in_tmp;  //!< Input samples for inverse FFT temporary storage
	wsgc_complex *_ifft_code_out_tmp; //!< Output samples for inverse FFT temporary storage
	unsigned int _prn_index;          //!< Current absolute PRN index
	std::vector<wsgc_float> _batch_sum_magnitudes;       //!< Sum of magnitudes at all iffti for the current batch and current prn
	std::vector<wsgc_float> _batch_noise_max_magnitude;  //!< Max of noise magnitude for the current batch
	unsigned int _buffer_multiplier;  //!< 1 for synced, 2 for non-synced (double buffer)


	/**
	 * Do the averaging when one batch has been processed
	 * \param message_correlation_records Reference to the message correlation records to which the results of the batch are appended
	 */
	void do_averaging(std::vector<CorrelationRecord>& message_correlation_records);

	/**
	 * Do the averaging when one batch has been processed - symbol synchronized version
	 * \param message_correlation_records Reference to the message correlation records to which the results of the batch are appended
	 */
	void do_averaging_synced(std::vector<CorrelationRecord>& message_correlation_records);
};

#endif /* __UNPILOTED_MESSAGE_CORRELATOR_H__ */
