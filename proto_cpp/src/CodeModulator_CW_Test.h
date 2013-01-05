/*
     Copyright 2012-2013 Edouard Griffiths <f4exb at free dot fr>
 
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
      
     CodeModulator for a continuous carrier. Basically does no perform any modulation.
     Used mainly for test signals in simulation in particular to test fading engines
      
*/
#ifndef __CODE_MODULATOR_CW_TEST_H__
#define __CODE_MODULATOR_CW_TEST_H__

#include "CodeModulator.h"

/**
 * \brief Special CodeModulator for CW signal. Basically does no perform any modulation.
 *  Used mainly for test signals in simulation in particular to test fading engines
 *
 * There is no difference btw fill_code_samples and modulate. Input signal is ignored. This is not a modulator actually
 * but is useful for simulation
 */
class CodeModulator_CW_Test : public CodeModulator
{
public:
	CodeModulator_CW_Test() {}
	~CodeModulator_CW_Test() {}
	virtual void fill_code_samples(wsgc_fftw_complex *fftw_code_in, std::vector<char>& code_bits);
	virtual void modulate(const wsgc_fftw_complex *in, wsgc_fftw_complex *out, std::vector<char>& code_bits);
};

#endif // __CODE_MODULATOR_CW_TEST_H__
