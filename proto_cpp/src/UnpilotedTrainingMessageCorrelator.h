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

     UnpilotedTrainingMessageCorrelator

     This flavour of correlator deals with training sequence PRNs without use of a pilot sequence

*/

#ifndef __UNPILOTED_TRAINING_MESSAGE_CORRELATOR_H__
#define __UNPILOTED_TRAINING_MESSAGE_CORRELATOR_H__

#include "WsgcTypes.h"
#include "TimeCorrelationAnalyzer.h"
#include <vector>

class TrainingCorrelationRecord;
/**
 * \brief Correlator engine to do correlation of PRN(s) in the training PRNs sequence without use of a pilot sequence
 */
class UnpilotedTrainingMessageCorrelator
{
public:
    /**
    * Correlator engine for training sequence PRNs
    * \param f_sampling Sampling frequency
    * \param f_chip Chip rate (frequency)
    * \param prn_length Length of a PRN sequence in number of chips
    * \param sequence_length Number of PRNs int the training sequence
    * \param averaging_length Number of PRNs used for sliding averaging
    * \param prn_list Reference to the vector of PRN numbers with which to make correlation
    * \param training_correlation_records Reference to the training correlation records
    */
	UnpilotedTrainingMessageCorrelator(
			wsgc_float f_sampling,
			wsgc_float f_chip,
			unsigned int prn_length,
            unsigned int sequence_length,
            unsigned int averaging_length,
			const std::vector<unsigned int>& prn_list,
			std::vector<TrainingCorrelationRecord>& training_correlation_records
			);
        
	virtual ~UnpilotedTrainingMessageCorrelator();

	/**
	 * Do the message correlation over the symbols window.
	 */
	virtual void execute() = 0;

	/**
	 * Append source samples for one PRN length to the buffer
     * \param samples Pointer to source samples
	 */
	virtual void set_samples(wsgc_complex *samples) = 0;

	/**
	 * Dump correlation records data to output stream
	 * \param os Output string stream
	 * \param mag_factor Magnitude correction factor
	 */
	void dump_correlation_records(std::ostringstream& os, wsgc_float mag_factor = 1.0);

	/**
	 * Dump message time shift analyzer results
	 */
	void dump_time_analyzer_results(std::ostringstream& os);

protected:
    wsgc_float _f_sampling; //!< Sampling frequency
    wsgc_float _f_chip; //!< Chip rate
    unsigned int _prn_length; //!< Length of a PRN sequence
    unsigned int _fft_N; //!< Length of FFT that is also the length of a PRN sequence in number of samples
    unsigned int _sequence_length; //!< Number of PRNs int the training sequence
    unsigned int _averaging_length; //!< Number of PRNs used for sliding averaging
    unsigned int _global_prn_index; //!< Index of PRN in all received samples
    const std::vector<unsigned int>& _prn_list; //!< Reference to the vector of PRN numbers with which to make correlation
    unsigned int _samples_length; //!< Number of stored source samples
    unsigned int _prn_in_avg_count; //!< Number of PRNs taking part in averaging
    unsigned int _prn_in_seq_count; //!< Number of PRNs taking part in sequence
	std::vector<TrainingCorrelationRecord>& _training_correlation_records; //!< Reference to the training correlation records
	TimeCorrelationAnalyzer<TrainingCorrelationRecord> _training_time_analyzer;
};

#endif /* __UNPILOTED_TRAINING_MESSAGE_CORRELATOR_H__ */
