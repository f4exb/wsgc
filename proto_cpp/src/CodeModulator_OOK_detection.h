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
      
     CodeModulator for detected baseband OOK signals
      
*/
#ifndef __CODE_MODULATOR_OOK_DETECTION_H__
#define __CODE_MODULATOR_OOK_DETECTION_H__

#include "CodeModulator.h"

/**
 * \brief Special CodeModulator for detected baseband OOK signals. Simulates an envelopped detected (real) OOK signal
 *
 * i.e. there is no difference btw fill_code_samples and modulate. Input signal is ignored. This is not a modulator actually
 * but is useful for simulation
 */
class CodeModulator_OOK_detection : public CodeModulator
{
public:
	CodeModulator_OOK_detection() {}
	virtual ~CodeModulator_OOK_detection() {}
	virtual void fill_code_samples(wsgc_fftw_complex *fftw_code_in, std::vector<char>& code_bits);
	virtual void modulate(const wsgc_fftw_complex *in, wsgc_fftw_complex *out, std::vector<char>& code_bits);
    virtual void fill_code_samples(wsgc_fftw_complex *fftw_code_in, std::vector<char>& code_bits, wsgc_float f_sampling, wsgc_float f_chip);
    virtual void modulate(const wsgc_fftw_complex *in, wsgc_fftw_complex *out, std::vector<char>& code_bits, wsgc_float f_sampling, wsgc_float f_chip);
};

#endif // __CODE_MODULATOR_OOK_DETECTION_H__
