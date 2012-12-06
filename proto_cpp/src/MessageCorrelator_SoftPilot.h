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

     MessageCorrelator

     Given the frequency and time displacement of correlation peak given by the
     Pilot Correlator it searches correlation for all PRNs of the message symbols
     to select the one that was sent. It uses straightforward time correlation.

*/

#ifndef __MESSAGE_CORRELATOR_H__
#define __MESSAGE_CORRELATOR_H__

#include "WsgcTypes.h"
#include "CorrelationRecord.h"
#include "ContinuousPhaseCarrier.h"
#include "MessageCorrelationMatrices.h"
#include <vector>

class GoldCodeGenerator ;
class LocalCodes;

/**
 * \brief Correlator engine to acquire and track message PRN(s) using the frequency
 * and time displacement of correlation peak given by the Pilot Correlator
 */
class MessageCorrelator_SoftPilot
{
public:
    /**
    * Correlator engine for message PRNs
    * \param gc_generator Gold Code generator being used
    * \param local_codes PRN signals local copy
    * \param f_sampling Sampling frequency
    * \param f_chip Chip rate (frequency)
    * \param prn_per_symbol Number of PRNs per symbol or averaging block
    * \param nb_batch_prns Number of PRNs processed in one batch ("PRN batch factor")
    */
	MessageCorrelator_SoftPilot(GoldCodeGenerator& gc_generator, LocalCodes& local_codes, wsgc_float f_sampling, wsgc_float f_chip, unsigned int prn_per_symbol, unsigned int nb_batch_prns);
	~MessageCorrelator_SoftPilot();

	/**
	 * Make a local copy of source FFT samples (one PRN length) applying frequency correction (multiplication)
	 * \param fft_src Pointer to the source FFT samples
     * \param delta_f Receiving frequency adjustment
	 */
	void copy_input_samples(wsgc_complex *fft_src, wsgc_float delta_f);

	/**
	 * Do the correlation on one PRN. 
     * \param t_shift PRN sequence time shift
     * \param global_prn_index PRN global sequence index count 
	 */
	void execute(unsigned int t_shift, unsigned int global_prn_index, CorrelationRecord& correlation_record);
    
protected:
    GoldCodeGenerator& _gc_generator; //!< Reference to the Gold Code generator being used
    LocalCodes& _local_codes; //!< Reference to the PRN signals local copy.
    ContinuousPhaseCarrier _local_oscillator; //!< Local oscillator for receiving frequency adjustment
    wsgc_float _f_sampling; //!< Sampling frequency
    wsgc_float _f_chip; //!< Chip rate
    unsigned int _nb_msg_prns; //!< Number of message PRNs to explore
    unsigned int _prn_per_symbol; //!< Number of PRNs per symbol that is also the number of PRNs per averaging block
    unsigned int _fft_N; //!< Size of FFT
    unsigned int _nb_batch_prns; //!< Number of PRNs in one batch of processing
    wsgc_float _delta_f; //!< Retain receiving frequency
    wsgc_complex *_src; //!< Source samples. Batch size.
    wsgc_complex *_corr_results; //!< Result of correlation of PRNs. Sparse matrix data part (PRNi x Pi)
    unsigned int *_corr_delta_t; //!< Delta t for the corresponding correlation result. Sparse matrix index part (index by delta t).
    wsgc_complex *_corr_avgsums; //!< Result of PRNs correlation averages. Sparse matrix PRNi x [(delta t, max)]times PRN per symbol
    std::vector<wsgc_float> _max_magnitudes; //!< Record maximum magnitude of all PRN correlations
    std::vector<unsigned int> _max_prn_indexes; //!< Record maximum (PRN) index of all PRN correlations
    MessageCorrelationMatrices _correlation_matrices; //!< Message correlation matrices manager
};

#endif /* __MESSAGE_CORRELATOR_H__ */
