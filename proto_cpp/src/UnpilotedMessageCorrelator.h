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

	 This is the super abstract class for Host and CUDA implementations

*/

#ifndef __UNPILOTED_MESSAGE_CORRELATOR_H__
#define __UNPILOTED_MESSAGE_CORRELATOR_H__

#include "WsgcTypes.h"
#include <vector>

class CorrelationRecord;

/**
 * \brief Correlator engine on all possible message PRNs (normally also including noise PRN)
 * to find a correlation peak representing the time delay of the start of sequence.
 *
 * Uses frequency domain correlation to check all possible delays at once which is
 * much more computationnaly efficient
 *
 * This is the super abstract class for Host and CUDA implementations
 */
class UnpilotedMessageCorrelator
{
public:
    /**
    * Correlator engine for message PRNs
    * \param f_sampling Sampling frequency
    * \param f_chip Chip rate (frequency)
    * \param prn_per_symbol Number of PRNs per symbol or averaging block
    * \param nb_batch_prns Number of PRNs per process batch
    */
	UnpilotedMessageCorrelator(
			wsgc_float f_sampling,
			wsgc_float f_chip,
			unsigned int prn_per_symbol,
			unsigned int nb_batch_prns
			) :
		_f_sampling(f_sampling),
		_f_chip(f_chip),
		_prn_per_symbol(prn_per_symbol),
		_nb_batch_prns(nb_batch_prns),
	    _batch_max_magnitudes(nb_batch_prns,0.0),
	    _batch_max_ifft_indexes(nb_batch_prns,0),
	    _batch_max_prn_indexes(nb_batch_prns,0),
	    _batch_complex_values_max(nb_batch_prns, (1.0, 0.0))
	{}

	virtual ~UnpilotedMessageCorrelator()
	{}

	/**
	 * Assign one PRN length of samples to be processed
	 * \param source_block Samples the length of one PRN to be processed
	 */
	virtual void set_source_block(wsgc_complex *source_block) = 0;

	/**
	 * Do the message correlation over the length of one analysis window.
	 * \param message_correlation_records Reference to the message correlation records that are the results of message correlation
     * \return True if a new averaging took plce hence new results are available
	 */
	virtual bool execute(std::vector<CorrelationRecord>& message_correlation_records) = 0;
    
    /**
     * Get the number of PRNs in one processing batch
     * \return The number of PRNs in one processing batch
     */
    unsigned int get_nb_batch_prns() const
    {
        return _nb_batch_prns;
    }
    
protected:
    wsgc_float _f_sampling; //!< Sampling frequency
    wsgc_float _f_chip; //!< Chip rate
    unsigned int _prn_per_symbol; //!< Number of PRNs per symbol that is also the number of PRNs per averaging block
    unsigned int _nb_batch_prns; //!< Number of PRNs per process batch

	std::vector<wsgc_float> _batch_max_magnitudes;       //!< Max magnitudes for the current batch
	std::vector<unsigned int> _batch_max_ifft_indexes;   //!< Corresponding IFFT indexes
	std::vector<unsigned int> _batch_max_prn_indexes;    //!< Corresponding PRN indexes
	std::vector<wsgc_complex> _batch_complex_values_max; //!< Max complex values for the current batch;

};

#endif /* __UNPILOTED_MESSAGE_CORRELATOR_H__ */
