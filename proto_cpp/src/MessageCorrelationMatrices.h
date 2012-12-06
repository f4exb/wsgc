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

     MessageCorrelationMatrices

	 This class encapsulates message correlation sparse matrices used with piloted resolution
	 The matrices are:
	 - Correlation results: 3 dimensional Delta-t x PRNi x Pi where Pi is the PRN sequence within symbol
	 	 This matrix has only one non null column (PRNi) per PRNi x Delta-t plane
	 	 It is stored as one PRNi x Pi matrix of correlation values and one PRNi vector of Delta-t values
	 - Averaging sum results: 2 dimensional PRNi x Delta-t
	 	 This matrix has only one non null column (PRNi) per PRNi x Delta-t plane
	 	 It is stored as one PRNi x Pi matrix of correlation values and one PRNi vector of Delta-t values
*/

#ifndef __MESSAGE_CORRELATION_MATRICES_H__
#define __MESSAGE_CORRELATION_MATRICES_H__

#include "WsgcTypes.h"
#include <vector>
#include <map>

/**
 * \brief This class encapsulates message correlation sparse matrices used with piloted resolution
 *   The matrices are:
	 - Correlation results: 3 dimensional Delta-t x PRNi x Pi where Pi is the PRN sequence within symbol
	 	 This matrix has only one non null column (PRNi) per PRNi x Delta-t plane
	 	 It is stored as one PRNi x Pi matrix of correlation values and one PRNi vector of Delta-t values
	 - Averaging sum results: 2 dimensional PRNi x Delta-t
	 	 This matrix has only one non null column (PRNi) per PRNi x Delta-t plane
	 	 It is stored as one PRNi x Pi matrix of correlation values and one PRNi vector of Delta-t values
 */
class MessageCorrelationMatrices
{
public:
	typedef struct correlation_tuple_s
	{
		unsigned int prni; //!< PRN index in PRN message alphabet
		unsigned int delta_t; //!< PRN start delay (Delta-t)
		wsgc_complex value; //!< Correlation value
	} CorrelationTuple_t;

    /**
    * MessageCorrelationMatrices mangement class
    * \param nb_message_prns Number of message PRNs to correlate plus one noise PRN (the last)
    * \param prn_per_symbol Number of PNRs per symbol
    * \param fft_N FFT size that is also the largest possible Delta-t
    */
	MessageCorrelationMatrices(unsigned int nb_message_prns, unsigned int prn_per_symbol, unsigned int fft_N);
	virtual ~MessageCorrelationMatrices();

	/**
	 * Add a correlation item
	 * \param prni PRN index in message alphabet
	 * \param pi PRN index in symbol
	 * \param delta_t PRN sequence delay
	 * \param value Correlation value
	 */
	void add_correlation_item(unsigned int prni, unsigned int pi, unsigned int delta_t, wsgc_complex value);

	/**
	 * Prepare for new PRNi vector
	 * \param pi PRN index in symbol
	 */
	void validate_prni_vector(unsigned int pi);

	/**
	 * Do an averaging process normally after each addition
	 */
	void process_averaging();

	/**
	 * Get the current correlation tuple with maximum value magnitude
	 * \param correlation_tuple Reference of the correlation tuple to be filled with result
	 * \param max_mag Maximum magnitude (of selected PRN)
     * \param avg_mag Average of maximum magnitudes of all PRNs
	 */
	void get_mag_max(CorrelationTuple_t& correlation_tuple, wsgc_float &max_mag, wsgc_float &avg_mag);

	/**
	 * Get the current correlation tuple with maximum noise value magnitude
	 * \param correlation_tuple Reference of the correlation tuple to be filled with result
	 * \param max_mag Maximum magnitude
	 */
	void get_noise_mag_max(CorrelationTuple_t& correlation_tuple, wsgc_float &max_mag);

	/**
	 * Get the noise average.
	 * \return noise average
	 */
	wsgc_float get_noise_avg();

protected:
	unsigned int _nb_message_prns; //!< Number of message PRNs to correlate plus one noise PRN (the last)
	unsigned int _prn_per_symbol; //!< Number of PNRs per symbol
	unsigned int _fft_N; //!< FFT size that is also the largest possible Delta-t
	wsgc_complex *_corr_values; //!< PRNi x Pi matrix of values
	unsigned int *_corr_delta_t; //!< Pi vector of Delta-t
	std::map<unsigned int, std::vector<wsgc_complex> > _avg_dict;
	wsgc_complex *_avg_values; //!< PRNi x Pi matrix of averaging values
	unsigned int *_avg_delta_t; //!< Pi vector of averaging Delta-t
	unsigned int _nb_averaging_sums; //!< Number of active (PRNi, Delta-t, value) tuples in _averaging_sums vector
	unsigned int _nb_corr_vectors; //!< Current number of correlation PRNi vectors stored. Increments at each storage until _prn_per_symbol is reached
	unsigned int _corr_item_index; //!< Current correlation item storage index
	unsigned int _corr_pi_at_zero_index; //!< PRN index in symbol cycle at start of correlation storage ring
};

#endif /* __MESSAGE_CORRELATION_MATRICES_H__ */
