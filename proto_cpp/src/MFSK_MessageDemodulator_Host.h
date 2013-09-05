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

     This is the Host implementation

*/
#ifndef __MFSK_MESSAGE_DEMODULATOR_HOST__
#define __MFSK_MESSAGE_DEMODULATOR_HOST__

#include "MFSK_MessageDemodulator.h"

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
class MFSK_MessageDemodulator_Host : public MFSK_MessageDemodulator
{
public:
	/**
	 * Build a new MFSK incoherent demodulator object
	 * \param fft_N Individual FFT size
	 * \param nb_fft_per_symbol Number of FFTs in one symbol
	 * \param zero_frequency_slot FFT frequency slot of ordinal zero symbol frequency
	 * \param nb_message_symbols Number of message symbols
	 * \param nb_service_symbols Number of service symbols (second is noise PRN)
	 */
	MFSK_MessageDemodulator_Host(
			unsigned int fft_N,
			unsigned int nb_fft_per_symbol,
			int zero_frequency_slot,
			unsigned int nb_message_symbols,
			unsigned int nb_service_symbols);

	virtual ~MFSK_MessageDemodulator_Host();

	/**
	 * Execute demodulation on one symbol length of samples. Time and frequency synchronization is supposed to have taken place
	 * Implementation (Host or CUDA) dependent
	 * \param symbol_samples Pointer to the symbol samples. Number of samples is assumed to be FFT size times the number of FFTs per symbol
     * \param relmat Reference of a RSSoft library reliability matrix.
	 */
#ifdef _RSSOFT
	virtual void execute(wsgc_complex *symbol_samples, rssoft::RS_ReliabilityMatrix& relmat);
#endif

    /**
     * Execute demodulation on one symbol length of samples. Time and frequency synchronization is supposed to have taken place
     * Implementation (Host or CUDA) dependent
     * \param symbol_samples Pointer to the symbol samples. Number of samples is assumed to be FFT size times the number of FFTs per symbol
     * \param relmat Pointer to a CCSoft library reliability matrix.
     */
#ifdef _CCSOFT
    virtual void execute(wsgc_complex *symbol_samples, ccsoft::CC_ReliabilityMatrix& relmat);
#endif

	virtual void execute(wsgc_complex *symbol_samples);

protected:
    wsgc_fftw_plan _fft_plan; //!< FFTW plan for forward FFT.
    wsgc_complex *_samples; //!< Source samples for one symbol length
    wsgc_complex *_src_fft; //!< Result of FFT of source samples for one symbol length
    wsgc_float *_magsum_s; //!< magnitude sum all FFTs in symbol by symbol index
    unsigned int _symbol_i; //!< index of symbol in message sequence
    unsigned int _fft_i; //!< index of FFT in symbol
    
    void clean_magsum();
    void cumulate_magsum(unsigned int fft_index);
    void estimate_magpeak();
    void calculate_magnitudes(wsgc_complex *symbol_samples);
};

#endif // __MFSK_MESSAGE_DEMODULATOR_HOST__
