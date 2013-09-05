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

     MFSK_MessageDemodulator

     Class to incoherently demodulate MFSK message. This is not using a correlation scheme
     and is there only for comparison to incoherent MFSK

     Abstract class to support Host or CUDA implementation

*/
#ifndef __MFSK_MESSAGE_DEMODULATOR__
#define __MFSK_MESSAGE_DEMODULATOR__

#include "WsgcTypes.h"
#include <vector>

class MFSK_MessageDemodulationRecord;

#ifdef _RSSOFT
namespace rssoft
{
    class RS_ReliabilityMatrix;
}
#endif

#ifdef _CCSOFT
namespace ccsoft
{
    class CC_ReliabilityMatrix;
}
#endif


/**
 * \brief Class to incoherently demodulate MFSK message. This is not using a correlation scheme
 *        and is there only for comparison to incoherent MFSK. Abstract class to support Host or CUDA implementation
 */
class MFSK_MessageDemodulator
{
public:
	/**
	 * Build a new MFSK incoherent demodulator object
	 * \param fft_N Individual FFT size
	 * \param nb_fft_per_symbol Number of FFTs in one symbol
	 * \param zero_fft_slot FFT slot of ordinal zero symbol
	 * \param nb_message_symbols Number of message symbols
	 * \param nb_service_symbols Number of service symbols (second is noise PRN)
	 */
	MFSK_MessageDemodulator(
			unsigned int _fft_N,
			unsigned int _nb_fft_per_symbol,
			int zero_ffty_slot,
			unsigned int nb_message_symbols,
			unsigned int nb_service_symbols);

	virtual ~MFSK_MessageDemodulator();

	/**
	 * Execute demodulation on one symbol length of samples. Time and frequency synchronization is supposed to have taken place
	 * Implementation (Host or CUDA) dependent
	 * \param symbol_samples Pointer to the symbol samples. Number of samples is assumed to be FFT size times the number of FFTs per symbol
     * \param relmat Reference of a RSSoft library reliability matrix.
	 */
#ifdef _RSSOFT
	virtual void execute(wsgc_complex *symbol_samples, rssoft::RS_ReliabilityMatrix& relmat) = 0;
#endif

    /**
     * Execute demodulation on one symbol length of samples. Time and frequency synchronization is supposed to have taken place
     * Implementation (Host or CUDA) dependent
     * \param symbol_samples Pointer to the symbol samples. Number of samples is assumed to be FFT size times the number of FFTs per symbol
     * \param relmat Pointer to a CCSoft library reliability matrix.
     */
#ifdef _CCSOFT
    virtual void execute(wsgc_complex *symbol_samples, ccsoft::CC_ReliabilityMatrix& relmat) = 0;
#endif

	/**
	 * Execute demodulation on one symbol length of samples. Time and frequency synchronization is supposed to have taken place
	 * Implementation (Host or CUDA) dependent
	 * \param symbol_samples Pointer to the symbol samples. Number of samples is assumed to be FFT size times the number of FFTs per symbol
	 */
    virtual void execute(wsgc_complex *symbol_samples) = 0;


	/**
	 * Dumps the demodulation records data to output stream
	 * \param os The output stream
	 * \param magnitude_factor Magnitudes are divided by this factor before display
	 */
	void dump_demodulation_records(std::ostringstream& os, wsgc_float magnitude_factor = 1.0) const;

	/**
	 * Get a reference to the message demodulation records
	 * \return Reference to the message demodulation records
	 */
	const std::vector<MFSK_MessageDemodulationRecord>& get_demodulation_records() const
	{
		return _demodulation_records;
	}

protected:
	unsigned int _fft_N; //!< Individual FFT size
	unsigned int _nb_fft_per_symbol; //!< Number of FFTs in one symbol
	int _zero_fft_slot; //!< FFT slot of ordinal zero symbol
	unsigned int _nb_message_symbols; //!< Number of message symbols
	unsigned int _nb_service_symbols; //!< Number of service symbols (second is noise PRN)
	std::vector<MFSK_MessageDemodulationRecord> _demodulation_records; //!< Demodulation results storage, one element per symbol

	int get_fft_slot(int symbol_ordinal) const;
	int get_symbol_ordinal(int fft_slot) const;
};

#endif // __MFSK_MESSAGE_DEMODULATOR__
