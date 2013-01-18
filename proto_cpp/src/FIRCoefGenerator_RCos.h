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

     FIRCoefGenerator_RCos

     FIR (real) coefficients generator for raised cosine filters

*/

#ifndef __FIR_COEF_GENERATOR_RCOS__
#define __FIR_COEF_GENERATOR_RCOS__

#include "FIRCoefGenerator.h"
#include <sstream>


/**
 * \brief FIR (real) coefficients generator for raised cosine filters
 */
class FIRCoefGenerator_RCos : public FIRCoefGenerator
{
public:
	/**
	 * Make a new set of tap coefficients
	 * \param f_sampling Sampling frequency
	 * \param f_cutoff Cutoff frequency. Must be twice the chip frequency for complex signals
	 * \param alpha Rolloff factor. Must be in [0.0,1.0] interval. Defaults to 1.0 otherwise.
	 * \param nb_taps Number of taps. Must be an odd number > 1. Defaults to next odd number if even.
	 */
	FIRCoefGenerator_RCos(
			wsgc_float f_sampling,
			wsgc_float f_cutoff,
			wsgc_float alpha,
			unsigned int nb_taps);

	virtual ~FIRCoefGenerator_RCos();

	/**
	 * Dump filter characteristics to output stream
	 */
	virtual void dump(std::ostringstream& os) const;


protected:
	wsgc_float _alpha;       //!< Rolloff factor in [0,1] interval

	static wsgc_float sinc(wsgc_float x);
	wsgc_float rcos(wsgc_float x);
};


#endif // __FIR_COEF_GENERATOR_RCOS__
