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
      
     Source FFT
      
     Multiplies input samples by each of the sub-frequency LOs
     Do the FFT of the result
     
*/
#ifndef __SOURCE_FFT__
#define __SOURCE_FFT__

#include "WsgcTypes.h"
#include <vector>

class CodeModulator;
class GoldCodeGenerator;

/**
 * \brief Local copy of codes to be used for frequency domain correlation
 *
   Creates a local copy of the complex conjugate of the FFT transform of all symbols gold codes
   This pre-calculates the code to be inserted in the final IFFT.
 *
 */
class SourceFFT
{
public:
	/**
	 * Constructor
	 * \param f_sampling Sampling frequency
	 * \param f_chip Chip rate
	 * \param fft_N FFT size
	 * \param freq_step_division Number of frequency sub steps
	 */
	SourceFFT(wsgc_float f_sampling,
			  wsgc_float f_chip,
			  unsigned int fft_N,
			  unsigned int freq_step_division);

	virtual ~SourceFFT();


protected:
	wsgc_float _f_sampling; //!< Sampling frequency
	wsgc_float _f_chip; //!< Chip rate
	unsigned int _fft_N; //!< Number of samples in one code length = FFT size
	unsigned int _freq_step_division; //!< Number of frequency sub-steps
};

#endif // __SOURCE_FFT__
