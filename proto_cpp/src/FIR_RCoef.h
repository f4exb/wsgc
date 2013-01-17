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

     FIR_RCoef

     Complex FIR filter with real coefficients

*/

#ifndef __FIR_RCOEF__
#define __FIR_RCOEF__

#include "WsgcTypes.h"
#include <vector>

/**
 * \brief Complex FIR filter with real tap coefficients
 */
class FIR_RCoef
{
public:
	/**
	 * Construct an new Complex FIR filter with real tap coefficients
	 * \param nb_taps Number of taps. Must be odd and > 1. Forced to be odd by adding one if even number is given.
	 */
	FIR_RCoef(std::vector<wsgc_float>& tap_coefs);
	virtual ~FIR_RCoef();

	/**
	 * Calculate a new sample
	 */
	wsgc_complex calc(wsgc_complex in);

protected:
	wsgc_complex *_delay_line;
	std::vector<wsgc_float>& _tap_coefs;
	unsigned int _last_index;
	unsigned int _nb_taps;
	static const wsgc_complex _c_zero;
};


#endif

