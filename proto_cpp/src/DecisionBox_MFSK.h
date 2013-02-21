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
      
     Decision box

     Final demodulated data analysis and message estimation for MFSK.
*/
#ifndef __DECISION_BOX_MFSK__
#define __DECISION_BOX_MFSK__

#include "DecisionRecord.h"
#include <vector>

class MFSK_MessageDemodulationRecord;

/**
 * \brief Demodulated data analysis and message estimation for MFSK.
 */
class DecisionBox_MFSK
{
public:

	typedef enum
	{
		decision_status_true_accept,
		decision_status_true_reject,
		decision_status_false_accept,
		decision_status_false_reject
	} decision_status_t;

	/**
	 * Constructor
	 * \param fft_size Size of the FFT
	 * \param nb_fft_per_symbol Number of FFTs in one symbol
	 */
	DecisionBox_MFSK(
			unsigned int fft_size,
			unsigned int nb_fft_per_symbol);

	virtual ~DecisionBox_MFSK();

	/**
	 * Returns a reference to the vector of decoded symbols
	 * \return Reference to the vector of decoded symbols
	 */
	const std::vector<int>& get_decoded_symbols() const
	{
		return _decoded_symbols;
	}

	/**
	 * Does a message symbols estimation
	 * \param demodulation_records Demodulation records from the MFSK demodulator
	 */
	void estimate_symbols(const std::vector<MFSK_MessageDemodulationRecord>& demodulation_records);

	/**
     * Print decision errors
     * \param os The output stream
     * \param original_symbols The symbols that were sent
	 * \param demodulation_records Demodulation records from the MFSK demodulator
     * \param no_trivial Do not print trivial results (true accept status and no_record decision)
     */
    void dump_decision_status(
    		std::ostringstream& os,
    		std::vector<unsigned int>& original_symbols,
    		const std::vector<MFSK_MessageDemodulationRecord>& demodulation_records,
    		bool no_trivial=true) const;

protected:
	unsigned int _fft_size; //!< Size of the FFT, this is also the number of samples in one PRN
	unsigned int _nb_fft_per_symbol; //!< Number of FFTs in one symbol
	std::vector<int> _decoded_symbols; //!< Symbols decoded, -1 for erasure
	static const wsgc_float peak_margin_threshold;  //!< selected correlation peak max difference with next / correlations ratio threshold

	void dump_decoding_status(std::ostringstream& os, decision_status_t decision_status) const;
};

#endif // __DECISION_BOX_MFSK__
