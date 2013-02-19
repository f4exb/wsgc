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
      
     CodeModulator for MFSK signals
      
*/
#ifndef __CODE_MODULATOR_MFSK_H__
#define __CODE_MODULATOR_MFSK_H__

#include "ContinuousPhaseCarrier.h"
#include <vector>

/**
 * \brief Code modulator for Multiple Frequency Shift Keying (MFSK). This is a special modulator
 * used for performance comparison with incoherent MFSK. It is not used in correlation.
 */
class CodeModulator_MFSK
{
public:
    /**
    * \brief Constructor
    * Makes a new MFSK code modulator object
    * \param f_sampling sampling frequency
    * \param zero_frequency Frequency of symbol zero (lowest frequency)
    * \param frequency_shift Symbol frequency shift (symbol bandwidth)
    * \param symbol_duration Symbol duration (symbol time)
    */
	CodeModulator_MFSK(
			wsgc_float f_sampling,
			wsgc_float zero_frequency,
			wsgc_float frequency_shift,
			wsgc_float symbol_duration);

	virtual ~CodeModulator_MFSK();
	/**
	 * Modulates the symbols into a waveform
	 * \param out pointer to the first element in the output waveform array. It must be allocated with at least of the size of the code_bits vector times the symbol length.
	 * \param symbols list of symbols
	 */
	void modulate(wsgc_fftw_complex *out, std::vector<unsigned int>& symbols);

	/**
	 * Returns the number of samples per symbol
	 * \return The number of samples per symbol
	 */
	unsigned int get_nb_symbol_samples() const
	{
		return _nb_symbol_samples;
	}

protected:
	ContinuousPhaseCarrier _local_oscillator;
	wsgc_float _zero_frequency;
	wsgc_float _frequency_shift;
	wsgc_float _symbol_duration;
	unsigned int _nb_symbol_samples;

};

#endif // __CODE_MODULATOR_MFSK_H__
