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

     This flavour of correlator deals with PRNs without use of a pilot sequence

*/

#ifndef __UNPILOTED_MULTIPLE_PRN_CORRELATOR_H__
#define __UNPILOTED_MULTIPLE_PRN_CORRELATOR_H__

#include "WsgcTypes.h"
#include "TimeCorrelationAnalyzer.h"
#include <vector>

class CorrelationRecord;

/**
 * \brief Correlator engine to do correlation of PRN(s) without use of a pilot sequence
 */
class UnpilotedMultiplePrnCorrelator
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
    * \param time_analysis_window_size Number of symbols used for time analysis.
    * \param correlation_records Reference to the correlation records
    */
	UnpilotedMultiplePrnCorrelator(
			wsgc_float f_sampling,
			wsgc_float f_chip,
			unsigned int prn_length,
            unsigned int prn_per_symbol,
			const std::vector<unsigned int>& prn_list,
			unsigned int symbol_window_size,
			unsigned int time_analysis_window_size,
			std::vector<CorrelationRecord>& correlation_records
			);
        
	virtual ~UnpilotedMultiplePrnCorrelator();

	/**
	 * Do the message correlation over the symbols window.
	 */
	virtual void execute() = 0;

	/**
	 * Append source samples for one PRN length to the buffer
     * \param samples Pointer to source samples
     * \return true if the buffer is complete and ready for analyze
	 */
	virtual bool set_samples(wsgc_complex *samples) = 0;

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
    unsigned int _prn_per_symbol; //!< Number of PRNs per message symbol
    unsigned int _global_prn_index; //!< Index of PRN in all received samples
    const std::vector<unsigned int>& _prn_list; //!< Reference to the vector of PRN numbers with which to make correlation
    unsigned int _symbol_window_size; //!<  Window of symbols in use for processing.
    unsigned int _time_analysis_window_size; //!<  Number of symbols used for time analysis.
    unsigned int _samples_length; //!< Number of stored source samples
    unsigned int _prns_length; //!< Number of stored PRNs
    std::vector<wsgc_float> _max_sy_mags; //!< Maximum magnitudes for each symbol
    std::vector<wsgc_float> _max_sy_mags_prni; //!< Maximum magnitudes sums (at each PRN) for each symbol
    std::vector<unsigned int> _max_sy_iffti; //!< Maximum ifft indexes at maximum magnitudes for each symbol
    std::vector<unsigned int> _max_sy_prni; //!< Maximum ifft indexes at maximum magnitudes for each symbol
	std::vector<CorrelationRecord>& _correlation_records; //!< Reference to the correlation records
	TimeCorrelationAnalyzer<CorrelationRecord> _message_time_analyzer;

    void init_results();
};

#endif /* __UNPILOTED_MULTIPLE_PRN_CORRELATOR_H__ */
